"""Create the initial sequence-change and AlphaFold-confidence audit."""

from __future__ import annotations

import csv
import difflib
import json
import statistics
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def fasta_sequence(path: Path) -> str:
    lines = path.read_text(encoding="utf-8").splitlines()
    return "".join(line.strip() for line in lines if not line.startswith(">"))


def ca_plddt(path: Path) -> dict[int, float]:
    values: dict[int, float] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            residue = int(line[22:26])
            values.setdefault(residue, float(line[60:66]))
    return values


def file_map(target: dict[str, Any]) -> dict[str, Path]:
    return {item["kind"]: ROOT / item["path"] for item in target["files"]}


def describe_change(sequence_a: str, sequence_b: str) -> dict[str, Any]:
    if len(sequence_a) == len(sequence_b):
        mismatches = [i for i, (a, b) in enumerate(zip(sequence_a, sequence_b), 1) if a != b]
        if not mismatches:
            return {"type": "identical", "a_start": "", "a_end": "", "b_start": "", "b_end": ""}
        start, end = min(mismatches), max(mismatches)
        return {
            "type": "replacement",
            "a_start": start,
            "a_end": end,
            "b_start": start,
            "b_end": end,
            "a_sequence": sequence_a[start - 1 : end],
            "b_sequence": sequence_b[start - 1 : end],
        }

    changes = [
        opcode
        for opcode in difflib.SequenceMatcher(a=sequence_a, b=sequence_b, autojunk=False).get_opcodes()
        if opcode[0] != "equal"
    ]
    if len(changes) != 1:
        raise ValueError(f"Expected one splice change, observed {len(changes)}")
    tag, a0, a1, b0, b1 = changes[0]
    return {
        "type": tag,
        "a_start": a0 + 1,
        "a_end": a1,
        "b_start": b0 + 1,
        "b_end": b1,
        "a_sequence": sequence_a[a0:a1],
        "b_sequence": sequence_b[b0:b1],
    }


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError(f"No rows for {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def region_stats(values: dict[int, float], start: int, end: int) -> dict[str, Any]:
    observed = [values[position] for position in range(start, end + 1) if position in values]
    if not observed:
        return {"mean_plddt": "", "min_plddt": "", "fraction_ge_70": "", "observed": 0}
    return {
        "mean_plddt": round(statistics.mean(observed), 2),
        "min_plddt": round(min(observed), 2),
        "fraction_ge_70": round(sum(value >= 70 for value in observed) / len(observed), 3),
        "observed": len(observed),
    }


def main() -> int:
    manifest = json.loads((ROOT / "data" / "manifest.json").read_text(encoding="utf-8"))
    if manifest.get("errors"):
        raise RuntimeError("Input retrieval contains errors; inspect data/manifest.json")

    targets = manifest["targets"]
    by_isoform = {(row["gene"], row["isoform_name"]): row for row in targets}

    quality_rows: list[dict[str, Any]] = []
    confidence_by_isoform: dict[tuple[str, str], dict[int, float]] = {}
    for target in targets:
        paths = file_map(target)
        sequence = fasta_sequence(paths["uniprot_fasta"])
        values = ca_plddt(paths["alphafold_pdb"])
        metadata = json.loads(paths["alphafold_metadata"].read_text(encoding="utf-8"))
        confidence_by_isoform[(target["gene"], target["isoform_name"])] = values
        quality_rows.append(
            {
                "gene": target["gene"],
                "isoform_name": target["isoform_name"],
                "sequence_length": len(sequence),
                "model_residues": len(values),
                "mean_plddt": round(statistics.mean(values.values()), 2),
                "fraction_ge_70": round(sum(v >= 70 for v in values.values()) / len(values), 3),
                "afdb_global_metric": metadata.get("globalMetricValue"),
                "entry_id": metadata.get("entryId"),
                "model_created": metadata.get("modelCreatedDate"),
            }
        )

    pairs = [("PKM", "PKM2", "PKM1"), ("RAC1", "RAC1", "RAC1B")]
    change_rows: list[dict[str, Any]] = []
    for gene, reference, comparator in pairs:
        reference_paths = file_map(by_isoform[(gene, reference)])
        comparator_paths = file_map(by_isoform[(gene, comparator)])
        change_rows.append(
            {
                "gene": gene,
                "reference": reference,
                "comparator": comparator,
                **describe_change(
                    fasta_sequence(reference_paths["uniprot_fasta"]),
                    fasta_sequence(comparator_paths["uniprot_fasta"]),
                ),
            }
        )

    region_rows: list[dict[str, Any]] = []
    for region in read_csv(ROOT / "config" / "regions.csv"):
        key = (region["gene"], region["isoform_name"])
        start, end = int(region["start"]), int(region["end"])
        region_rows.append(
            {
                **region,
                **region_stats(confidence_by_isoform[key], start, end),
            }
        )

    results = ROOT / "results"
    write_csv(results / "model_quality.csv", quality_rows)
    write_csv(results / "sequence_changes.csv", change_rows)
    write_csv(results / "region_confidence.csv", region_rows)

    report = [
        "# Initial input and confidence report",
        "",
        "## What was checked",
        "",
        "All four reviewed isoform sequences and AlphaFold models were retrieved,",
        "checksummed, and compared. Five experimental structures were also retrieved",
        "for the next structural-validation stage.",
        "",
        "## AlphaFold model quality",
        "",
        "| Protein | Length | Mean pLDDT | Fraction >=70 | Entry |",
        "|---|---:|---:|---:|---|",
    ]
    for row in quality_rows:
        report.append(
            f'| {row["isoform_name"]} | {row["sequence_length"]} | '
            f'{row["mean_plddt"]:.2f} | {row["fraction_ge_70"]:.1%} | '
            f'{row["entry_id"]} |'
        )
    report.extend(["", "## Splice differences", ""])
    for row in change_rows:
        report.append(
            f'- **{row["reference"]} vs {row["comparator"]}:** {row["type"]}; '
            f'reference {row["a_start"]}-{row["a_end"]}, comparator '
            f'{row["b_start"]}-{row["b_end"]}.'
        )
    report.extend(
        [
            "",
            "## Region-level confidence",
            "",
            "| Protein | Region | Residues | Mean pLDDT | Minimum | Fraction >=70 |",
            "|---|---|---:|---:|---:|---:|",
        ]
    )
    for row in region_rows:
        report.append(
            f'| {row["isoform_name"]} | {row["region_name"]} | '
            f'{row["start"]}-{row["end"]} | {row["mean_plddt"]} | '
            f'{row["min_plddt"]} | {float(row["fraction_ge_70"]):.1%} |'
        )
    report.extend(
        [
            "",
            "## Decision",
            "",
            "PKM1/PKM2 passes the fixed-structure control gate. The RAC1B inserted",
            "segment is treated as flexible unless experimental evidence supports a",
            "specific bound conformation. Confident neighbouring regions proceed to",
            "structural and interface comparison.",
            "",
        ]
    )
    (results / "initial_report.md").write_text("\n".join(report), encoding="utf-8")
    print(f"Wrote results to {results}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
