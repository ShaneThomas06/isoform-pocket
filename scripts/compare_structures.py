"""Superimpose AlphaFold and experimental structures using conserved residues."""

from __future__ import annotations

import csv
import difflib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M",
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def parse_pdb_chains(path: Path) -> dict[str, list[dict[str, Any]]]:
    chains: dict[str, list[dict[str, Any]]] = {}
    seen: set[tuple[str, str, str]] = set()
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.startswith("ATOM") or line[12:16].strip() != "CA":
            continue
        chain = line[21].strip() or "_"
        residue_id = (chain, line[22:26].strip(), line[26].strip())
        if residue_id in seen:
            continue
        seen.add(residue_id)
        residue_name = line[17:20].strip()
        if residue_name not in AA3_TO_1:
            continue
        chains.setdefault(chain, []).append(
            {
                "chain": chain,
                "number": int(line[22:26]),
                "insertion_code": line[26].strip(),
                "aa": AA3_TO_1[residue_name],
                "coord": np.array(
                    [float(line[30:38]), float(line[38:46]), float(line[46:54])],
                    dtype=float,
                ),
            }
        )
    return chains


def longest_chain(path: Path) -> list[dict[str, Any]]:
    chains = parse_pdb_chains(path)
    if not chains:
        raise ValueError(f"No protein CA atoms in {path}")
    return max(chains.values(), key=len)


def sequence_map(a: list[dict[str, Any]], b: list[dict[str, Any]]) -> list[tuple[int, int]]:
    sequence_a = "".join(row["aa"] for row in a)
    sequence_b = "".join(row["aa"] for row in b)
    mapping: list[tuple[int, int]] = []
    matcher = difflib.SequenceMatcher(a=sequence_a, b=sequence_b, autojunk=False)
    for tag, a0, a1, b0, b1 in matcher.get_opcodes():
        if tag == "equal":
            mapping.extend(zip(range(a0, a1), range(b0, b1)))
    return mapping


def kabsch(mobile: np.ndarray, reference: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    mobile_center = mobile.mean(axis=0)
    reference_center = reference.mean(axis=0)
    covariance = (mobile - mobile_center).T @ (reference - reference_center)
    u, _, vt = np.linalg.svd(covariance)
    if np.linalg.det(u @ vt) < 0:
        u[:, -1] *= -1
    rotation = u @ vt
    translation = reference_center - mobile_center @ rotation
    return rotation, translation


def rmsd(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.sum((a - b) ** 2, axis=1))))


def superposition_metrics(
    reference: list[dict[str, Any]],
    mobile: list[dict[str, Any]],
    mapping: list[tuple[int, int]] | None = None,
) -> tuple[dict[str, Any], np.ndarray, np.ndarray, list[tuple[int, int]]]:
    mapping = sequence_map(reference, mobile) if mapping is None else mapping
    if len(mapping) < 3:
        raise ValueError("Fewer than three conserved residues could be mapped")
    reference_coords = np.array([reference[i]["coord"] for i, _ in mapping])
    mobile_coords = np.array([mobile[j]["coord"] for _, j in mapping])
    rotation, translation = kabsch(mobile_coords, reference_coords)
    aligned_mobile = mobile_coords @ rotation + translation
    distances = np.linalg.norm(reference_coords - aligned_mobile, axis=1)
    metrics = {
        "mapped_residues": len(mapping),
        "rmsd_angstrom": round(rmsd(reference_coords, aligned_mobile), 3),
        "median_displacement": round(float(np.median(distances)), 3),
        "p90_displacement": round(float(np.percentile(distances, 90)), 3),
        "max_displacement": round(float(np.max(distances)), 3),
    }
    return metrics, rotation, translation, mapping


def manifest_paths(target: dict[str, Any]) -> dict[str, Path]:
    return {item["kind"]: ROOT / item["path"] for item in target["files"]}


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    manifest = json.loads((ROOT / "data" / "manifest.json").read_text(encoding="utf-8"))
    targets = {(row["gene"], row["isoform_name"]): row for row in manifest["targets"]}
    experimental = {
        row["pdb_id"]: ROOT / row["files"][0]["path"]
        for row in manifest["experimental_structures"]
    }

    comparisons: list[dict[str, Any]] = []
    displacement_rows: list[dict[str, Any]] = []

    af_pairs = [
        ("PKM", "PKM2", "PKM1"),
        ("RAC1", "RAC1", "RAC1B"),
    ]
    for gene, reference_name, mobile_name in af_pairs:
        reference = longest_chain(manifest_paths(targets[(gene, reference_name)])["alphafold_pdb"])
        mobile = longest_chain(manifest_paths(targets[(gene, mobile_name)])["alphafold_pdb"])
        if gene == "PKM":
            # The two PKM isoforms have equal length and direct positional
            # equivalence. Exclude the entire alternative-exon region from the
            # fit, then measure it without letting it define the reference frame.
            alignment_mapping = [
                (index, index) for index in range(len(reference)) if not 388 <= index <= 432
            ]
            displacement_mapping = [(index, index) for index in range(len(reference))]
        else:
            alignment_mapping = sequence_map(reference, mobile)
            displacement_mapping = alignment_mapping
        metrics, rotation, translation, _ = superposition_metrics(
            reference, mobile, alignment_mapping
        )
        comparisons.append(
            {
                "gene": gene,
                "reference": f"AF:{reference_name}",
                "mobile": f"AF:{mobile_name}",
                "comparison_type": "alphafold_isoform_pair",
                **metrics,
            }
        )
        for reference_index, mobile_index in displacement_mapping:
            ref_coord = reference[reference_index]["coord"]
            mobile_coord = mobile[mobile_index]["coord"] @ rotation + translation
            displacement_rows.append(
                {
                    "gene": gene,
                    "reference_isoform": reference_name,
                    "mobile_isoform": mobile_name,
                    "reference_residue": reference[reference_index]["number"],
                    "mobile_residue": mobile[mobile_index]["number"],
                    "amino_acid": reference[reference_index]["aa"],
                    "displacement_angstrom": round(float(np.linalg.norm(ref_coord - mobile_coord)), 3),
                }
            )

    validation_pairs = [
        ("PKM", "PKM1", "3SRF"),
        ("PKM", "PKM2", "3U2Z"),
        ("RAC1", "RAC1", "8S1N"),
        ("RAC1", "RAC1B", "1RYF"),
        ("RAC1", "RAC1B", "1RYH"),
    ]
    for gene, isoform, pdb_id in validation_pairs:
        predicted = longest_chain(manifest_paths(targets[(gene, isoform)])["alphafold_pdb"])
        observed = longest_chain(experimental[pdb_id])
        metrics, _, _, _ = superposition_metrics(observed, predicted)
        comparisons.append(
            {
                "gene": gene,
                "reference": f"PDB:{pdb_id}",
                "mobile": f"AF:{isoform}",
                "comparison_type": "experimental_validation",
                **metrics,
            }
        )

    results = ROOT / "results"
    write_csv(results / "structure_comparisons.csv", comparisons)
    write_csv(results / "per_residue_displacement.csv", displacement_rows)

    def region_summary(gene: str, ref_start: int, ref_end: int) -> tuple[float, float, int]:
        values = [
            float(row["displacement_angstrom"])
            for row in displacement_rows
            if row["gene"] == gene and ref_start <= int(row["reference_residue"]) <= ref_end
        ]
        return round(statistics_mean(values), 3), round(max(values), 3), len(values)

    pkm_local = region_summary("PKM", 389, 433)
    rac_upstream = region_summary("RAC1", 66, 75)
    rac_downstream = region_summary("RAC1", 76, 85)

    report = [
        "# Structural comparison report",
        "",
        "Structures were aligned only on sequence-identical residues using a Kabsch",
        "least-squares superposition. The splice-altered residues therefore did not",
        "determine the alignment frame.",
        "",
        "## Global comparisons",
        "",
        "| Gene | Reference | Compared model | Type | Mapped residues | RMSD (A) | P90 displacement (A) |",
        "|---|---|---|---|---:|---:|---:|",
    ]
    for row in comparisons:
        report.append(
            f'| {row["gene"]} | {row["reference"]} | {row["mobile"]} | '
            f'{row["comparison_type"]} | {row["mapped_residues"]} | '
            f'{row["rmsd_angstrom"]} | {row["p90_displacement"]} |'
        )
    report.extend(
        [
            "",
            "## Isoform-local displacement after global superposition",
            "",
            f"- PKM reference residues 389-433: mean {pkm_local[0]} A; maximum {pkm_local[1]} A; n={pkm_local[2]} positionally equivalent residue pairs.",
            f"- RAC1 residues 66-75, immediately upstream of the insertion: mean {rac_upstream[0]} A; maximum {rac_upstream[1]} A; n={rac_upstream[2]}.",
            f"- RAC1 residues 76-85, mapped to RAC1B residues 95-104 downstream of the insertion: mean {rac_downstream[0]} A; maximum {rac_downstream[1]} A; n={rac_downstream[2]}.",
            "",
            "## Interpretation boundary",
            "",
            "A displacement is a geometric observation, not proof of a druggable site.",
            "The next stage must test whether any displacement changes a reproducible",
            "cavity or known interaction surface and whether experimental structures",
            "support the same conclusion.",
            "",
        ]
    )
    (results / "structural_report.md").write_text("\n".join(report), encoding="utf-8")
    print(f"Wrote structural comparisons to {results}")
    return 0


def statistics_mean(values: list[float]) -> float:
    if not values:
        return math.nan
    return sum(values) / len(values)


if __name__ == "__main__":
    raise SystemExit(main())
