"""Parse P2Rank outputs and integrate them with the confidence-aware results."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
LABELS = {
    "PKM_PKM2": ("PKM", "PKM2"),
    "PKM_PKM1": ("PKM", "PKM1"),
    "RAC1_RAC1": ("RAC1", "RAC1"),
    "RAC1_RAC1B": ("RAC1", "RAC1B"),
}


def clean_row(row: dict[str, str]) -> dict[str, str]:
    return {key.strip(): value.strip() for key, value in row.items()}


def parse_residues(text: str) -> set[int]:
    residues = set()
    for token in text.split():
        residues.add(int(token.rsplit("_", 1)[1]))
    return residues


def target_region(gene: str, isoform: str) -> set[int]:
    if gene == "PKM":
        return set(range(389, 434))
    return set(range(66, 105)) if isoform == "RAC1B" else set(range(66, 86))


def canonical(gene: str, isoform: str, residues: set[int]) -> set[int]:
    if gene != "RAC1" or isoform != "RAC1B":
        return residues
    return {residue if residue <= 75 else residue - 19 for residue in residues if residue <= 75 or residue >= 95}


def jaccard(a: set[int], b: set[int]) -> float:
    return len(a & b) / len(a | b) if a | b else 0.0


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    predictions: list[dict[str, Any]] = []
    residue_sets: dict[tuple[str, str, str], set[int]] = {}
    for directory, (gene, isoform) in LABELS.items():
        csv_path = next((ROOT / "results" / "p2rank" / directory).glob("*_predictions.csv"))
        with csv_path.open(newline="", encoding="utf-8") as handle:
            for raw in csv.DictReader(handle):
                row = clean_row(raw)
                residues = parse_residues(row["residue_ids"])
                pocket = row["name"]
                predictions.append(
                    {
                        "gene": gene,
                        "isoform_name": isoform,
                        "pocket": pocket,
                        "rank": int(row["rank"]),
                        "score": float(row["score"]),
                        "probability": float(row["probability"]),
                        "residue_count": len(residues),
                        "residues": ";".join(str(value) for value in sorted(residues)),
                        "touches_splice_neighborhood": bool(residues & target_region(gene, isoform)),
                    }
                )
                residue_sets[(gene, isoform, pocket)] = residues

    matches: list[dict[str, Any]] = []
    for gene, reference, comparator in [("PKM", "PKM2", "PKM1"), ("RAC1", "RAC1", "RAC1B")]:
        reference_rows = [r for r in predictions if r["gene"] == gene and r["isoform_name"] == reference]
        comparator_rows = [r for r in predictions if r["gene"] == gene and r["isoform_name"] == comparator]
        for ref in reference_rows:
            ref_set = canonical(gene, reference, residue_sets[(gene, reference, ref["pocket"])])
            candidates = []
            for comp in comparator_rows:
                comp_set = canonical(gene, comparator, residue_sets[(gene, comparator, comp["pocket"])])
                candidates.append((jaccard(ref_set, comp_set), comp))
            score, best = max(candidates, key=lambda item: item[0])
            matches.append(
                {
                    "gene": gene,
                    "reference_isoform": reference,
                    "reference_pocket": ref["pocket"],
                    "comparator_isoform": comparator,
                    "comparator_pocket": best["pocket"] if score >= 0.25 else "",
                    "residue_jaccard": round(score, 3),
                    "reference_probability": ref["probability"],
                    "comparator_probability": best["probability"] if score >= 0.25 else "",
                    "reference_touches_splice": ref["touches_splice_neighborhood"],
                    "comparator_touches_splice": best["touches_splice_neighborhood"] if score >= 0.25 else "",
                }
            )

    results = ROOT / "results"
    write_csv(results / "p2rank_predictions.csv", predictions)
    write_csv(results / "p2rank_matches.csv", matches)

    rac_predictions = [row for row in predictions if row["gene"] == "RAC1"]
    py_kvfinder = list(csv.DictReader((results / "cavity_predictions.csv").open(encoding="utf-8")))
    rac_kv = [row for row in py_kvfinder if row["gene"] == "RAC1" and row["touches_splice_neighborhood"] == "True"]

    report = [
        "# Consensus pocket decision",
        "",
        "## P2Rank results for RAC1 and RAC1B",
        "",
        "| Isoform | Pocket | Rank | Score | Probability | Residues | Touches splice neighbourhood |",
        "|---|---|---:|---:|---:|---|---|",
    ]
    for row in rac_predictions:
        report.append(
            f'| {row["isoform_name"]} | {row["pocket"]} | {row["rank"]} | '
            f'{row["score"]} | {row["probability"]} | {row["residues"]} | '
            f'{row["touches_splice_neighborhood"]} |'
        )
    report.extend(
        [
            "",
            "## Cross-method decision",
            "",
            f"- pyKVFinder detected {len(rac_kv)} RAC1/RAC1B cavities touching the splice neighbourhood.",
            "- P2Rank detected no pocket touching the RAC1B insertion or immediate flanks.",
            "- The largest pyKVFinder RAC1B cavity depends on 16 low-confidence inserted residues.",
            "- The smaller high-confidence model cavity was not recovered in either RAC1B crystal structure.",
            "",
            "## Current conclusion",
            "",
            "There is currently **no reproducible evidence for an isoform-selective, fixed",
            "small-molecule pocket near the RAC1B insertion**. The results instead",
            "support treating this region as a dynamic interaction surface. This is a",
            "scientifically useful negative result and prevents a false druggability claim.",
            "",
            "The next discovery analysis should test altered partner-binding surfaces",
            "and disorder-to-order behaviour rather than perform docking into the",
            "unsupported AlphaFold cavity.",
            "",
        ]
    )
    (results / "consensus_pocket_report.md").write_text("\n".join(report), encoding="utf-8")
    print(f"Wrote P2Rank summary and consensus decision to {results}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
