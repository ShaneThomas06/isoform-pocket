"""Analyze matched RAC1 and RAC1B AlphaFold Server jobs across seeds."""

from __future__ import annotations

import csv
import json
import os
import re
import statistics
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLBACKEND", "Agg")

ROOT = Path(__file__).resolve().parents[1]
AF_ROOT = ROOT / "data" / "alphafold_server"
RESULTS = ROOT / "results"
FIGURES = RESULTS / "figures"
SEEDS = (202, 303, 404)
ISOFORMS = ("RAC1B", "RAC1")
HOTSPOT = (395, 430, 433, 434, 437, 474, 477)
CTNND1_OFFSET = 349
CONTACT_THRESHOLD = 0.30


def model_number(path: Path) -> int:
    match = re.search(r"_(\d+)\.json$", path.name)
    if not match:
        raise ValueError(f"Cannot parse model number from {path.name}")
    return int(match.group(1))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"No rows generated for {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def job_directory(isoform: str, seed: int) -> Path:
    return AF_ROOT / f"{isoform}_p120ARM350-824_seed{seed}" / "extracted"


def required_files(directory: Path, pattern: str) -> list[Path]:
    files = sorted(directory.glob(pattern), key=model_number)
    if len(files) != 5:
        raise FileNotFoundError(
            f"Expected five files matching {pattern} in {directory}, found {len(files)}"
        )
    return files


def mean_matrix(matrices: list[list[list[float]]]) -> list[list[float]]:
    size = len(matrices[0])
    return [
        [statistics.mean(matrix[i][j] for matrix in matrices) for j in range(size)]
        for i in range(size)
    ]


def analyze_job(isoform: str, seed: int) -> tuple[list[dict[str, Any]], dict[str, Any], list[dict[str, Any]]]:
    directory = job_directory(isoform, seed)
    summary_files = required_files(directory, "*_summary_confidences_*.json")
    full_files = required_files(directory, "*_full_data_*.json")

    model_rows: list[dict[str, Any]] = []
    summaries: list[dict[str, Any]] = []
    for path in summary_files:
        row = json.loads(path.read_text(encoding="utf-8"))
        summaries.append(row)
        model_rows.append(
            {
                "isoform": isoform,
                "seed": seed,
                "model": model_number(path),
                "iptm": row["iptm"],
                "ptm": row["ptm"],
                "ranking_score": row["ranking_score"],
                "pair_pae_min": row["chain_pair_pae_min"][0][1],
            }
        )

    full = [json.loads(path.read_text(encoding="utf-8")) for path in full_files]
    reference = full[0]
    chain_ids = reference["token_chain_ids"]
    residue_ids = reference["token_res_ids"]
    if any(item["token_chain_ids"] != chain_ids for item in full[1:]):
        raise ValueError(f"Token-chain mapping differs among models in {directory}")
    if any(item["token_res_ids"] != residue_ids for item in full[1:]):
        raise ValueError(f"Token-residue mapping differs among models in {directory}")

    contacts = mean_matrix([item["contact_probs"] for item in full])
    rac_indices = [i for i, chain in enumerate(chain_ids) if chain == "A"]
    p120_indices = [i for i, chain in enumerate(chain_ids) if chain == "B"]
    if isoform == "RAC1B":
        source_indices = [i for i in rac_indices if 76 <= residue_ids[i] <= 94]
        source_label = "RAC1B insertion 76-94"
    else:
        source_indices = rac_indices
        source_label = "RAC1 all residues"

    hotspot_rows: list[dict[str, Any]] = []
    for canonical_position in HOTSPOT:
        local_position = canonical_position - CTNND1_OFFSET
        matches = [
            i
            for i in p120_indices
            if residue_ids[i] == local_position
        ]
        if len(matches) != 1:
            raise ValueError(
                f"Expected one p120 token for CTNND1 {canonical_position}, found {len(matches)}"
            )
        p120_index = matches[0]
        value = max(contacts[rac_index][p120_index] for rac_index in source_indices)
        hotspot_rows.append(
            {
                "isoform": isoform,
                "seed": seed,
                "ctnnd1_position": canonical_position,
                "max_mean_contact_probability": round(value, 4),
                "source_residues": source_label,
            }
        )

    pair_count = sum(
        contacts[rac_index][p120_index] >= CONTACT_THRESHOLD
        for rac_index in rac_indices
        for p120_index in p120_indices
    )
    summary = {
        "isoform": isoform,
        "seed": seed,
        "mean_iptm": round(statistics.mean(row["iptm"] for row in summaries), 4),
        "mean_ptm": round(statistics.mean(row["ptm"] for row in summaries), 4),
        "pairs_ge_0_30": pair_count,
        "mean_hotspot_signal": round(
            statistics.mean(
                float(row["max_mean_contact_probability"])
                for row in hotspot_rows
            ),
            4,
        ),
    }
    return model_rows, summary, hotspot_rows


def write_report(summary_rows: list[dict[str, Any]], hotspot_rows: list[dict[str, Any]]) -> None:
    by_key = {
        (row["isoform"], int(row["seed"]), int(row["ctnnd1_position"])): float(
            row["max_mean_contact_probability"]
        )
        for row in hotspot_rows
    }
    lines = [
        "# Multi-seed RAC1B/RAC1-p120 ARM-domain analysis",
        "",
        "## Conclusion",
        "",
        "The RAC1B-specific contact pattern observed at seed 202 does not reproduce",
        "at seeds 303 and 404. RAC1 also contacts the same p120 surface. The data do",
        "not support an isoform-exclusive interface.",
        "",
        "Every seed-level mean ipTM is below 0.23. No predicted docking orientation",
        "is interpreted as a resolved complex.",
        "",
        "## Seed-level summary",
        "",
        "| Isoform | Seed | Mean ipTM | Mean pTM | All-chain pairs >= 0.30 | Mean hotspot signal |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in summary_rows:
        lines.append(
            f"| {row['isoform']} | {row['seed']} | {float(row['mean_iptm']):.3f} | "
            f"{float(row['mean_ptm']):.3f} | {row['pairs_ge_0_30']} | "
            f"{float(row['mean_hotspot_signal']):.3f} |"
        )
    lines.extend(
        [
            "",
            "## Across-seed candidate-surface averages",
            "",
            "| CTNND1 residue | RAC1B insertion | RAC1 whole protein | Difference |",
            "|---:|---:|---:|---:|",
        ]
    )
    for position in HOTSPOT:
        rac1b = statistics.mean(by_key[("RAC1B", seed, position)] for seed in SEEDS)
        rac1 = statistics.mean(by_key[("RAC1", seed, position)] for seed in SEEDS)
        lines.append(f"| {position} | {rac1b:.3f} | {rac1:.3f} | {rac1b-rac1:+.3f} |")
    lines.extend(
        [
            "",
            "## Claim boundary",
            "",
            "The candidate surface is seed-sensitive and shared by both isoforms.",
            "Small average differences at individual residues remain descriptive.",
            "They are not treated as binding determinants without orthogonal data.",
            "",
        ]
    )
    (RESULTS / "af3_multiseed_report.md").write_text("\n".join(lines), encoding="utf-8")


def make_figure(summary_rows: list[dict[str, Any]], hotspot_rows: list[dict[str, Any]]) -> None:
    import matplotlib.pyplot as plt
    import numpy as np

    colors = {"RAC1B": "#ea6a0a", "RAC1": "#0f6b8d"}
    seed_colors = {202: "#7042c1", 303: "#2a9d8f", 404: "#d99a0b"}
    by_summary = {(row["isoform"], int(row["seed"])): row for row in summary_rows}
    by_hotspot = {
        (row["isoform"], int(row["seed"]), int(row["ctnnd1_position"])): float(
            row["max_mean_contact_probability"]
        )
        for row in hotspot_rows
    }

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.titlesize": 15,
            "axes.labelsize": 12,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
        }
    )
    fig, axes = plt.subplots(2, 2, figsize=(16, 11), facecolor="white")
    fig.subplots_adjust(left=0.08, right=0.96, bottom=0.08, top=0.82, wspace=0.30, hspace=0.45)
    fig.suptitle("Multi-seed reproducibility: RAC1B and RAC1 with the p120 ARM domain", y=0.96, fontsize=23, fontweight="bold")
    fig.text(0.5, 0.915, "Seeds 202, 303, and 404; five diffusion models per job", ha="center", fontsize=14, color="#475569")

    ax = axes[0, 0]
    for isoform in ISOFORMS:
        ax.plot(SEEDS, [float(by_summary[(isoform, seed)]["mean_iptm"]) for seed in SEEDS], marker="o", linewidth=2.5, markersize=7, color=colors[isoform], label=isoform)
    ax.axhspan(0, 0.3, color="#fee2e2", alpha=0.65)
    ax.set_ylim(0, 0.5)
    ax.set_xticks(SEEDS)
    ax.set_xlabel("Seed")
    ax.set_ylabel("Mean ipTM")
    ax.set_title("A  Global complex confidence remains low", loc="left", pad=10)
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.2)

    ax = axes[0, 1]
    rows = [(isoform, seed) for isoform in ISOFORMS for seed in SEEDS]
    matrix = np.array([[by_hotspot[(isoform, seed, position)] for position in HOTSPOT] for isoform, seed in rows])
    image = ax.imshow(matrix, cmap="YlOrRd", vmin=0, vmax=0.65, aspect="auto")
    ax.set_xticks(range(len(HOTSPOT)), labels=HOTSPOT)
    ax.set_yticks(range(len(rows)), labels=[f"{isoform} {seed}" for isoform, seed in rows])
    ax.set_xlabel("Canonical CTNND1 residue")
    ax.set_title("B  Candidate-surface signal varies by seed", loc="left", pad=10)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            ax.text(j, i, f"{matrix[i,j]:.2f}", ha="center", va="center", fontsize=9, color="white" if matrix[i,j] >= 0.36 else "#2b2118", fontweight="bold")
    colorbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.03)
    colorbar.set_label("Maximum mean contact probability")

    ax = axes[1, 0]
    for isoform in ISOFORMS:
        ax.plot(SEEDS, [float(by_summary[(isoform, seed)]["mean_hotspot_signal"]) for seed in SEEDS], marker="o", linewidth=2.5, markersize=7, color=colors[isoform], label=isoform)
    ax.set_xticks(SEEDS)
    ax.set_ylim(0, 0.55)
    ax.set_xlabel("Seed")
    ax.set_ylabel("Mean signal across seven positions")
    ax.set_title("C  Seed-202 specificity does not replicate", loc="left", pad=10)
    ax.grid(axis="y", alpha=0.2)

    ax = axes[1, 1]
    for seed in SEEDS:
        delta = [by_hotspot[("RAC1B", seed, position)] - by_hotspot[("RAC1", seed, position)] for position in HOTSPOT]
        ax.plot(HOTSPOT, delta, marker="o", linewidth=2.2, markersize=6, color=seed_colors[seed], label=f"Seed {seed}")
    ax.axhline(0, color="#334155", linewidth=1.2)
    ax.set_xlabel("Canonical CTNND1 residue")
    ax.set_ylabel("RAC1B minus RAC1 contact probability")
    ax.set_title("D  Isoform differences change with the seed", loc="left", pad=10)
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.2)

    FIGURES.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURES / "figure7_af3_multiseed_reproducibility.png", dpi=220, facecolor="white")
    plt.close(fig)


def main() -> int:
    all_models: list[dict[str, Any]] = []
    summaries: list[dict[str, Any]] = []
    hotspots: list[dict[str, Any]] = []
    for isoform in ISOFORMS:
        for seed in SEEDS:
            model_rows, summary, hotspot_rows = analyze_job(isoform, seed)
            all_models.extend(model_rows)
            summaries.append(summary)
            hotspots.extend(hotspot_rows)

    RESULTS.mkdir(parents=True, exist_ok=True)
    write_csv(RESULTS / "af3_multiseed_model_metrics.csv", all_models)
    write_csv(RESULTS / "af3_multiseed_summary.csv", summaries)
    write_csv(RESULTS / "af3_multiseed_hotspot_values.csv", hotspots)
    write_report(summaries, hotspots)
    make_figure(summaries, hotspots)
    print(f"Analyzed {len(all_models)} AlphaFold models across {len(summaries)} matched jobs.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
