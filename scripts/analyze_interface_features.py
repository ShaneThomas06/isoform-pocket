"""Analyze sequence chemistry and model exposure around the RAC1B insertion."""

from __future__ import annotations

import csv
import json
import math
import statistics
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np

try:
    from scripts.analyze_inputs import fasta_sequence
    from scripts.compare_structures import longest_chain
except ModuleNotFoundError:  # Support direct execution from the repository root.
    from analyze_inputs import fasta_sequence
    from compare_structures import longest_chain


ROOT = Path(__file__).resolve().parents[1]
POSITIVE = set("KR")
NEGATIVE = set("DE")
HYDROPHOBIC = set("AVILMFWY")
AROMATIC = set("FWY")
FLEXIBLE = set("GPS")
HYDROPATHY = {
    "A": 1.8,
    "R": -4.5,
    "N": -3.5,
    "D": -3.5,
    "C": 2.5,
    "Q": -3.5,
    "E": -3.5,
    "G": -0.4,
    "H": -3.2,
    "I": 4.5,
    "L": 3.8,
    "K": -3.9,
    "M": 1.9,
    "F": 2.8,
    "P": -1.6,
    "S": -0.8,
    "T": -0.7,
    "W": -0.9,
    "Y": -1.3,
    "V": 4.2,
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"No rows for {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def sequence_features(sequence: str) -> dict[str, Any]:
    length = len(sequence)
    counts = Counter(sequence)
    entropy = -sum(
        (count / length) * math.log2(count / length) for count in counts.values()
    )
    return {
        "length": length,
        "net_charge_proxy": sum(counts[aa] for aa in POSITIVE)
        - sum(counts[aa] for aa in NEGATIVE),
        "positive_fraction": round(sum(counts[aa] for aa in POSITIVE) / length, 3),
        "negative_fraction": round(sum(counts[aa] for aa in NEGATIVE) / length, 3),
        "hydrophobic_fraction": round(
            sum(counts[aa] for aa in HYDROPHOBIC) / length, 3
        ),
        "aromatic_fraction": round(sum(counts[aa] for aa in AROMATIC) / length, 3),
        "flexible_fraction": round(sum(counts[aa] for aa in FLEXIBLE) / length, 3),
        "mean_hydropathy": round(
            statistics.mean(HYDROPATHY[aa] for aa in sequence), 3
        ),
        "sequence_entropy": round(entropy, 3),
    }


def exposure_by_residue(pdb_path: Path) -> dict[int, dict[str, float]]:
    """Calculate a C-alpha neighbour exposure proxy, not a formal SASA."""

    residues = longest_chain(pdb_path)
    coords = np.array([row["coord"] for row in residues])
    output: dict[int, dict[str, float]] = {}
    for index, row in enumerate(residues):
        distances = np.linalg.norm(coords - coords[index], axis=1)
        neighbours = sum(
            distance <= 10.0 and abs(other - index) > 2
            for other, distance in enumerate(distances)
        )
        output[row["number"]] = {
            "ca_neighbours_10a": int(neighbours),
            "exposure_proxy": round(1.0 / (1.0 + neighbours), 4),
        }
    return output


def plddt_by_residue(pdb_path: Path) -> dict[int, float]:
    output: dict[int, float] = {}
    for line in pdb_path.read_text(encoding="utf-8").splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            output.setdefault(int(line[22:26]), float(line[60:66]))
    return output


def region_label(isoform: str, position: int) -> str:
    if isoform == "RAC1B":
        if 66 <= position <= 75:
            return "upstream_flank"
        if 76 <= position <= 94:
            return "insertion"
        if 95 <= position <= 104:
            return "downstream_flank"
    return "matched_neighbourhood"


def chemistry_class(aa: str) -> str:
    if aa in POSITIVE:
        return "positive"
    if aa in NEGATIVE:
        return "negative"
    if aa in AROMATIC:
        return "aromatic"
    if aa in HYDROPHOBIC:
        return "hydrophobic"
    return "polar_or_flexible"


def target_paths(
    targets: dict[tuple[str, str], dict[str, Any]], isoform: str
) -> dict[str, Path]:
    row = targets[("RAC1", isoform)]
    return {item["kind"]: ROOT / item["path"] for item in row["files"]}


def candidate_rows(residue_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    insertion = [
        row
        for row in residue_rows
        if row["isoform"] == "RAC1B" and row["region"] == "insertion"
    ]
    candidates: list[dict[str, Any]] = []
    for row in insertion:
        aa = row["amino_acid"]
        rationale: list[str] = []
        priority = 0.0
        if aa in NEGATIVE:
            priority += 3.0
            rationale.append("tests recognition of the p120 polybasic loop")
        if aa in POSITIVE:
            priority += 2.0
            rationale.append("tests the insertion's basic cluster")
        if aa in AROMATIC:
            priority += 2.5
            rationale.append("tests a possible transient aromatic anchor")
        if aa == "P":
            priority += 1.5
            rationale.append("tests insertion conformational preference")
        if not rationale:
            continue
        priority += row["exposure_proxy"]
        candidates.append(
            {
                "isoform": "RAC1B",
                "position": row["position"],
                "amino_acid": aa,
                "suggested_substitution": "A",
                "priority_score": round(priority, 3),
                "rationale": "; ".join(rationale),
                "confidence_warning": (
                    "Low insertion pLDDT: hypothesis-generating, not a fixed interface."
                ),
            }
        )
    return sorted(
        candidates,
        key=lambda row: (-float(row["priority_score"]), int(row["position"])),
    )


def main() -> int:
    manifest = json.loads((ROOT / "data" / "manifest.json").read_text(encoding="utf-8"))
    targets = {
        (row["gene"], row["isoform_name"]): row for row in manifest["targets"]
    }

    region_rows: list[dict[str, Any]] = []
    for row in read_csv(ROOT / "config" / "interface_regions.csv"):
        region_rows.append(
            {
                "protein": row["protein"],
                "isoform": row["isoform"],
                "start": row["start"],
                "end": row["end"],
                "sequence": row["sequence"],
                **sequence_features(row["sequence"]),
            }
        )

    residue_rows: list[dict[str, Any]] = []
    for isoform, start, end in (("RAC1", 66, 85), ("RAC1B", 66, 104)):
        paths = target_paths(targets, isoform)
        sequence = fasta_sequence(paths["uniprot_fasta"])
        confidence = plddt_by_residue(paths["alphafold_pdb"])
        exposure = exposure_by_residue(paths["alphafold_pdb"])
        for position in range(start, end + 1):
            aa = sequence[position - 1]
            residue_rows.append(
                {
                    "isoform": isoform,
                    "position": position,
                    "amino_acid": aa,
                    "region": region_label(isoform, position),
                    "chemistry": chemistry_class(aa),
                    "charge_proxy": (
                        1 if aa in POSITIVE else -1 if aa in NEGATIVE else 0
                    ),
                    "hydropathy": HYDROPATHY[aa],
                    "plddt": confidence[position],
                    **exposure[position],
                    "exposure_interpretation": (
                        "exploratory_low_confidence"
                        if confidence[position] < 70
                        else "usable_model_proxy"
                    ),
                }
            )

    candidates = candidate_rows(residue_rows)
    results = ROOT / "results"
    write_csv(results / "interface_region_features.csv", region_rows)
    write_csv(results / "interface_residue_features.csv", residue_rows)
    write_csv(results / "interface_candidate_residues.csv", candidates)

    rac1b, p120 = region_rows
    report = [
        "# RAC1B-p120 interface feature report",
        "",
        "## Numbering and evidence",
        "",
        "The binding experiments used short p120-catenin isoform 4. This analysis",
        "retains canonical human CTNND1/O60716 coordinates for the shared ARM-domain",
        "sequence (607-644). RAC1B uses P63000-2 coordinates. No complex structure",
        "is assumed.",
        "",
        "## Sequence-level comparison",
        "",
        "| Region | Sequence | Net charge | Positive | Negative | Hydrophobic | Mean hydropathy |",
        "|---|---|---:|---:|---:|---:|---:|",
        (
            f'| RAC1B 76-94 | `{rac1b["sequence"]}` | {rac1b["net_charge_proxy"]} | '
            f'{float(rac1b["positive_fraction"]):.1%} | '
            f'{float(rac1b["negative_fraction"]):.1%} | '
            f'{float(rac1b["hydrophobic_fraction"]):.1%} | '
            f'{rac1b["mean_hydropathy"]} |'
        ),
        (
            f'| p120 607-644 | `{p120["sequence"]}` | {p120["net_charge_proxy"]} | '
            f'{float(p120["positive_fraction"]):.1%} | '
            f'{float(p120["negative_fraction"]):.1%} | '
            f'{float(p120["hydrophobic_fraction"]):.1%} | '
            f'{p120["mean_hydropathy"]} |'
        ),
        "",
        "Both segments have mixed charge and flexible residues. Because the p120",
        "segment is net positive, RAC1B acidic residues are direct experimental",
        "candidates, but composition alone cannot specify a binding pose.",
        "",
        "## First-pass mutagenesis candidates",
        "",
        "| Rank | RAC1B residue | Suggested test | Reason |",
        "|---:|---|---|---|",
    ]
    for rank, row in enumerate(candidates[:5], 1):
        report.append(
            f'| {rank} | {row["amino_acid"]}{row["position"]} | '
            f'{row["amino_acid"]}{row["position"]}A | {row["rationale"]} |'
        )
    report.extend(
        [
            "",
            "## Interpretation",
            "",
            "These candidates test electrostatic, aromatic, and conformational",
            "features of the insertion. The ranking is not a prediction of",
            "experimental effect. Because insertion pLDDT is low, AlphaFold exposure",
            "values are exploratory metadata and do not define a fixed interface.",
            "",
        ]
    )
    (results / "interface_feature_report.md").write_text(
        "\n".join(report), encoding="utf-8"
    )
    print("Wrote interface feature tables and report.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
