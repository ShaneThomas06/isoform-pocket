"""Check AlphaFold cavities against monomeric experimental structures."""

from __future__ import annotations

import csv
import json
import tempfile
from pathlib import Path
from typing import Any

import pyKVFinder

from compare_structures import longest_chain, sequence_map


ROOT = Path(__file__).resolve().parents[1]


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def residue_set(text: str) -> set[int]:
    return {int(value) for value in text.split(";") if value}


def jaccard(a: set[int], b: set[int]) -> float:
    union = a | b
    return len(a & b) / len(union) if union else 0.0


def protein_only_chain(source: Path, destination: Path, chain: str) -> None:
    lines = []
    for line in source.read_text(encoding="utf-8").splitlines():
        if line.startswith("ATOM") and (line[21].strip() or "_") == chain:
            alternate = line[16]
            if alternate not in {" ", "A"}:
                continue
            if alternate == "A":
                line = line[:16] + " " + line[17:]
            lines.append(line)
    lines.extend(["TER", "END"])
    destination.write_text("\n".join(lines) + "\n", encoding="utf-8")


def file_map(target: dict[str, Any]) -> dict[str, Path]:
    return {item["kind"]: ROOT / item["path"] for item in target["files"]}


def main() -> int:
    manifest = json.loads((ROOT / "data" / "manifest.json").read_text(encoding="utf-8"))
    targets = {(row["gene"], row["isoform_name"]): row for row in manifest["targets"]}
    predicted_rows = read_csv(ROOT / "results" / "cavity_predictions.csv")

    validation_rows: list[dict[str, Any]] = []
    experimental_rows: list[dict[str, Any]] = []

    with tempfile.TemporaryDirectory(prefix="isoform-pocket-") as temporary:
        temporary_dir = Path(temporary)
        for record in manifest["experimental_structures"]:
            gene, isoform, pdb_id = record["gene"], record["isoform_name"], record["pdb_id"]
            source = ROOT / record["files"][0]["path"]
            experimental_chain = longest_chain(source)
            chain_id = experimental_chain[0]["chain"]
            cleaned = temporary_dir / f"{pdb_id}_{chain_id}.pdb"
            protein_only_chain(source, cleaned, chain_id)

            predicted_path = file_map(targets[(gene, isoform)])["alphafold_pdb"]
            predicted_chain = longest_chain(predicted_path)
            mapping = sequence_map(experimental_chain, predicted_chain)
            experimental_to_af = {
                experimental_chain[exp_index]["number"]: predicted_chain[af_index]["number"]
                for exp_index, af_index in mapping
            }

            print(f"Detecting experimental cavities in {pdb_id} chain {chain_id}", flush=True)
            result = pyKVFinder.run_workflow(
                str(cleaned),
                step=0.6,
                probe_in=1.4,
                probe_out=4.0,
                removal_distance=2.4,
                volume_cutoff=5.0,
                include_depth=True,
                include_hydropathy=True,
                nthreads=4,
            )

            predicted_for_isoform = [
                row
                for row in predicted_rows
                if row["gene"] == gene and row["isoform_name"] == isoform
            ]
            for cavity_id in sorted(result.volume):
                experimental_numbers = {int(row[0]) for row in result.residues[cavity_id]}
                mapped_numbers = {
                    experimental_to_af[number]
                    for number in experimental_numbers
                    if number in experimental_to_af
                }
                candidates = [
                    (jaccard(mapped_numbers, residue_set(row["residues"])), row)
                    for row in predicted_for_isoform
                ]
                score, best = max(candidates, key=lambda item: item[0])
                experimental_rows.append(
                    {
                        "gene": gene,
                        "isoform_name": isoform,
                        "pdb_id": pdb_id,
                        "chain": chain_id,
                        "experimental_cavity": cavity_id,
                        "volume_a3": round(float(result.volume[cavity_id]), 3),
                        "max_depth_a": round(float(result.max_depth[cavity_id]), 3),
                        "mapped_residues": ";".join(str(value) for value in sorted(mapped_numbers)),
                        "mapped_residue_count": len(mapped_numbers),
                        "best_af_cavity": best["cavity_id"],
                        "residue_jaccard": round(score, 3),
                    }
                )
                if score >= 0.25:
                    validation_rows.append(
                        {
                            "gene": gene,
                            "isoform_name": isoform,
                            "af_cavity": best["cavity_id"],
                            "af_volume_a3": best["volume_a3"],
                            "pdb_id": pdb_id,
                            "experimental_cavity": cavity_id,
                            "experimental_volume_a3": round(float(result.volume[cavity_id]), 3),
                            "residue_jaccard": round(score, 3),
                            "af_touches_splice": best["touches_splice_neighborhood"],
                            "af_inserted_residue_count": best["inserted_residue_count"],
                        }
                    )

    results = ROOT / "results"
    write_csv(results / "experimental_cavities.csv", experimental_rows)
    write_csv(results / "cavity_validation.csv", validation_rows)

    rac_relevant = [
        row
        for row in validation_rows
        if row["gene"] == "RAC1" and row["af_touches_splice"] == "True"
    ]
    report = [
        "# Experimental cavity validation",
        "",
        "Each experimental PDB was reduced to its longest protein chain and all",
        "HETATM records were removed before cavity detection. Experimental residue",
        "numbers were sequence-mapped back to the corresponding AlphaFold model.",
        "This validates monomeric cavities only, not oligomer-interface sites.",
        "",
        "## RAC1/RAC1B splice-neighbourhood matches",
        "",
        "| Isoform | AF cavity | AF volume (A3) | PDB | Experimental cavity | Experimental volume (A3) | Residue Jaccard | Inserted residues in AF cavity |",
        "|---|---|---:|---|---|---:|---:|---:|",
    ]
    for row in rac_relevant:
        report.append(
            f'| {row["isoform_name"]} | {row["af_cavity"]} | {row["af_volume_a3"]} | '
            f'{row["pdb_id"]} | {row["experimental_cavity"]} | '
            f'{row["experimental_volume_a3"]} | {row["residue_jaccard"]} | '
            f'{row["af_inserted_residue_count"]} |'
        )
    report.extend(
        [
            "",
            "## Decision rule",
            "",
            "A RAC1B cavity is retained only if it has meaningful experimental",
            "residue overlap and does not depend primarily on the unresolved insertion.",
            "The result still requires an independent pocket detector before being",
            "described as a supported isoform-selective candidate.",
            "",
        ]
    )
    (results / "experimental_validation_report.md").write_text(
        "\n".join(report), encoding="utf-8"
    )
    print(f"Wrote experimental cavity validation to {results}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
