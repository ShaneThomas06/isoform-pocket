"""Detect cavities in AlphaFold models and compare isoform-adjacent sites."""

from __future__ import annotations

import csv
import json
import statistics
from pathlib import Path
from typing import Any

import pyKVFinder


ROOT = Path(__file__).resolve().parents[1]


def file_map(target: dict[str, Any]) -> dict[str, Path]:
    return {item["kind"]: ROOT / item["path"] for item in target["files"]}


def ca_plddt(path: Path) -> dict[int, float]:
    values: dict[int, float] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            values.setdefault(int(line[22:26]), float(line[60:66]))
    return values


def interface_residue_numbers(rows: list[list[str]]) -> set[int]:
    return {int(row[0]) for row in rows}


def canonical_residues(gene: str, isoform: str, residues: set[int]) -> set[int]:
    if gene != "RAC1" or isoform != "RAC1B":
        return set(residues)
    normalized: set[int] = set()
    for residue in residues:
        if residue <= 75:
            normalized.add(residue)
        elif residue >= 95:
            normalized.add(residue - 19)
    return normalized


def target_region(gene: str, isoform: str) -> set[int]:
    if gene == "PKM":
        return set(range(389, 434))
    if isoform == "RAC1B":
        return set(range(66, 105))
    return set(range(66, 86))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"No rows generated for {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def jaccard(a: set[int], b: set[int]) -> float:
    union = a | b
    return len(a & b) / len(union) if union else 0.0


def main() -> int:
    manifest = json.loads((ROOT / "data" / "manifest.json").read_text(encoding="utf-8"))
    cavity_rows: list[dict[str, Any]] = []
    residue_sets: dict[tuple[str, str, str], set[int]] = {}

    for target in manifest["targets"]:
        gene, isoform = target["gene"], target["isoform_name"]
        pdb_path = file_map(target)["alphafold_pdb"]
        print(f"Detecting cavities in {gene}/{isoform}", flush=True)
        result = pyKVFinder.run_workflow(
            str(pdb_path),
            step=0.6,
            probe_in=1.4,
            probe_out=4.0,
            removal_distance=2.4,
            volume_cutoff=5.0,
            include_depth=True,
            include_hydropathy=True,
            nthreads=4,
        )
        confidence = ca_plddt(pdb_path)
        region = target_region(gene, isoform)
        for cavity_id in sorted(result.volume):
            residues = interface_residue_numbers(result.residues[cavity_id])
            local_confidence = [confidence[number] for number in residues if number in confidence]
            inserted_residues = (
                len(residues & set(range(76, 95)))
                if gene == "RAC1" and isoform == "RAC1B"
                else 0
            )
            row = {
                "gene": gene,
                "isoform_name": isoform,
                "cavity_id": cavity_id,
                "volume_a3": round(float(result.volume[cavity_id]), 3),
                "area_a2": round(float(result.area[cavity_id]), 3),
                "mean_depth_a": round(float(result.avg_depth[cavity_id]), 3),
                "max_depth_a": round(float(result.max_depth[cavity_id]), 3),
                "mean_hydropathy": round(float(result.avg_hydropathy[cavity_id]), 3),
                "residue_count": len(residues),
                "residues": ";".join(str(value) for value in sorted(residues)),
                "touches_splice_neighborhood": bool(residues & region),
                "inserted_residue_count": inserted_residues,
                "mean_residue_plddt": round(statistics.mean(local_confidence), 2),
                "min_residue_plddt": round(min(local_confidence), 2),
                "fraction_residues_ge_70": round(
                    sum(value >= 70 for value in local_confidence) / len(local_confidence), 3
                ),
            }
            cavity_rows.append(row)
            residue_sets[(gene, isoform, cavity_id)] = residues

    pair_definitions = [("PKM", "PKM2", "PKM1"), ("RAC1", "RAC1", "RAC1B")]
    match_rows: list[dict[str, Any]] = []
    for gene, reference_isoform, comparator_isoform in pair_definitions:
        reference_rows = [
            row for row in cavity_rows if row["gene"] == gene and row["isoform_name"] == reference_isoform
        ]
        comparator_rows = [
            row for row in cavity_rows if row["gene"] == gene and row["isoform_name"] == comparator_isoform
        ]
        used_comparator: set[str] = set()
        for reference in reference_rows:
            ref_set = canonical_residues(
                gene,
                reference_isoform,
                residue_sets[(gene, reference_isoform, reference["cavity_id"])],
            )
            candidates: list[tuple[float, dict[str, Any]]] = []
            for comparator in comparator_rows:
                comp_set = canonical_residues(
                    gene,
                    comparator_isoform,
                    residue_sets[(gene, comparator_isoform, comparator["cavity_id"])],
                )
                candidates.append((jaccard(ref_set, comp_set), comparator))
            score, best = max(candidates, key=lambda item: item[0])
            matched = score >= 0.25
            if matched:
                used_comparator.add(best["cavity_id"])
            match_rows.append(
                {
                    "gene": gene,
                    "reference_isoform": reference_isoform,
                    "reference_cavity": reference["cavity_id"],
                    "comparator_isoform": comparator_isoform,
                    "comparator_cavity": best["cavity_id"] if matched else "",
                    "residue_jaccard": round(score, 3),
                    "reference_volume_a3": reference["volume_a3"],
                    "comparator_volume_a3": best["volume_a3"] if matched else "",
                    "volume_change_a3": (
                        round(float(best["volume_a3"]) - float(reference["volume_a3"]), 3)
                        if matched
                        else ""
                    ),
                    "reference_touches_splice": reference["touches_splice_neighborhood"],
                    "comparator_touches_splice": (
                        best["touches_splice_neighborhood"] if matched else ""
                    ),
                    "preliminary_match": matched,
                }
            )
        for comparator in comparator_rows:
            if comparator["cavity_id"] not in used_comparator:
                match_rows.append(
                    {
                        "gene": gene,
                        "reference_isoform": reference_isoform,
                        "reference_cavity": "",
                        "comparator_isoform": comparator_isoform,
                        "comparator_cavity": comparator["cavity_id"],
                        "residue_jaccard": 0.0,
                        "reference_volume_a3": "",
                        "comparator_volume_a3": comparator["volume_a3"],
                        "volume_change_a3": "",
                        "reference_touches_splice": "",
                        "comparator_touches_splice": comparator["touches_splice_neighborhood"],
                        "preliminary_match": False,
                    }
                )

    results_dir = ROOT / "results"
    write_csv(results_dir / "cavity_predictions.csv", cavity_rows)
    write_csv(results_dir / "cavity_matches.csv", match_rows)

    relevant = [row for row in cavity_rows if row["touches_splice_neighborhood"]]
    relevant.sort(key=lambda row: (row["gene"], -float(row["volume_a3"])))
    report = [
        "# Preliminary cavity report",
        "",
        "Cavities were detected with pyKVFinder 0.9.2 using the preregistered",
        "default geometry. This is the first detector only; independent replication",
        "is required before a strong pocket claim.",
        "",
        "## Cavities touching splice-site neighbourhoods",
        "",
        "| Gene | Isoform | Cavity | Volume (A3) | Max depth (A) | Mean pLDDT | Minimum pLDDT | Inserted residues |",
        "|---|---|---|---:|---:|---:|---:|---:|",
    ]
    for row in relevant:
        report.append(
            f'| {row["gene"]} | {row["isoform_name"]} | {row["cavity_id"]} | '
            f'{row["volume_a3"]} | {row["max_depth_a"]} | '
            f'{row["mean_residue_plddt"]} | {row["min_residue_plddt"]} | '
            f'{row["inserted_residue_count"]} |'
        )
    report.extend(
        [
            "",
            "## Interpretation",
            "",
            "Cavities involving low-confidence RAC1B insertion residues are flagged",
            "as exploratory. Priority is given to sites supported by confident flanking",
            "residues, an equivalent experimental cavity, and a second detector.",
            "",
        ]
    )
    (results_dir / "cavity_report.md").write_text("\n".join(report), encoding="utf-8")
    print(f"Wrote cavity analysis to {results_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
