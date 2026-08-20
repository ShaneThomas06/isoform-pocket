"""Generate publication-style figures for the first analysis milestone."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
from matplotlib.patches import Patch


ROOT = Path(__file__).resolve().parents[1]
COLORS = {"PKM2": "#225ea8", "PKM1": "#41b6c4", "RAC1": "#238b45", "RAC1B": "#d95f0e"}


def file_map(target: dict[str, Any]) -> dict[str, Path]:
    return {item["kind"]: ROOT / item["path"] for item in target["files"]}


def ca_plddt(path: Path) -> dict[int, float]:
    values: dict[int, float] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            values.setdefault(int(line[22:26]), float(line[60:66]))
    return values


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def main() -> int:
    figure_dir = ROOT / "results" / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((ROOT / "data" / "manifest.json").read_text(encoding="utf-8"))

    fig, axes = plt.subplots(2, 2, figsize=(12, 6.8), constrained_layout=True)
    for axis, target in zip(axes.flat, manifest["targets"]):
        isoform = target["isoform_name"]
        values = ca_plddt(file_map(target)["alphafold_pdb"])
        axis.plot(list(values), list(values.values()), color=COLORS[isoform], linewidth=1.25)
        axis.axhline(70, color="#555555", linestyle="--", linewidth=0.9, label="pLDDT 70")
        if target["gene"] == "PKM":
            axis.axvspan(389, 433, color="#756bb1", alpha=0.16, label="splice-altered region")
        elif isoform == "RAC1B":
            axis.axvspan(76, 94, color="#de2d26", alpha=0.18, label="inserted segment")
        else:
            axis.axvspan(66, 85, color="#756bb1", alpha=0.13, label="comparison region")
        axis.set(title=f'{target["gene"]}: {isoform}', xlabel="Residue", ylabel="pLDDT", ylim=(0, 102))
        axis.spines[["top", "right"]].set_visible(False)
        axis.legend(frameon=False, fontsize=8, loc="lower left")
    fig.suptitle("AlphaFold confidence is high except within the RAC1B insertion", fontsize=14)
    fig.savefig(figure_dir / "figure1_confidence_profiles.png", dpi=200)
    plt.close(fig)

    displacement = read_csv(ROOT / "results" / "per_residue_displacement.csv")
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.2), constrained_layout=True)
    settings = [
        ("PKM", "PKM2 vs PKM1", (389, 433), "#225ea8"),
        ("RAC1", "RAC1 vs RAC1B", (66, 85), "#d95f0e"),
    ]
    for axis, (gene, title, region, color) in zip(axes, settings):
        rows = [row for row in displacement if row["gene"] == gene]
        x = [int(row["reference_residue"]) for row in rows]
        y = [float(row["displacement_angstrom"]) for row in rows]
        axis.plot(x, y, color=color, linewidth=1.1)
        axis.scatter(x, y, color=color, s=7, alpha=0.65)
        axis.axvspan(region[0], region[1], color="#756bb1", alpha=0.14)
        axis.set(title=title, xlabel="Reference residue", ylabel="CA displacement after alignment (A)")
        axis.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Isoform differences concentrate at splice boundaries", fontsize=14)
    fig.savefig(figure_dir / "figure2_structural_displacement.png", dpi=200)
    plt.close(fig)

    p2rank = [row for row in read_csv(ROOT / "results" / "p2rank_predictions.csv") if row["gene"] == "RAC1"]
    labels = [f'{row["isoform_name"]}\n{row["pocket"]}' for row in p2rank]
    probabilities = [float(row["probability"]) for row in p2rank]
    colors = [COLORS[row["isoform_name"]] for row in p2rank]
    fig, axis = plt.subplots(figsize=(7.5, 4.2), constrained_layout=True)
    bars = axis.bar(labels, probabilities, color=colors, width=0.62)
    axis.axhline(0.5, color="#555555", linestyle="--", linewidth=0.9, label="0.5 probability")
    axis.set(ylabel="P2Rank pocket probability", ylim=(0, 0.62))
    axis.set_title("P2Rank finds no splice-neighbourhood pocket")
    axis.spines[["top", "right"]].set_visible(False)
    axis.legend(frameon=False)
    for bar, probability in zip(bars, probabilities):
        axis.text(bar.get_x() + bar.get_width() / 2, probability + 0.018, f"{probability:.3f}", ha="center")
    fig.savefig(figure_dir / "figure3_p2rank_probabilities.png", dpi=200)
    plt.close(fig)

    interface = [
        row
        for row in read_csv(ROOT / "results" / "interface_residue_features.csv")
        if row["isoform"] == "RAC1B" and row["region"] == "insertion"
    ]
    candidates = read_csv(ROOT / "results" / "interface_candidate_residues.csv")[:5]
    candidate_positions = {int(row["position"]) for row in candidates}
    chemistry_colors = {
        "positive": "#2b8cbe",
        "negative": "#de2d26",
        "aromatic": "#756bb1",
        "hydrophobic": "#31a354",
        "polar_or_flexible": "#969696",
    }
    positions = [int(row["position"]) for row in interface]
    labels = [f'{row["amino_acid"]}{row["position"]}' for row in interface]
    hydropathy = [float(row["hydropathy"]) for row in interface]
    confidence = [float(row["plddt"]) for row in interface]
    bar_colors = [chemistry_colors[row["chemistry"]] for row in interface]

    fig, axes = plt.subplots(
        2,
        1,
        figsize=(11, 6.2),
        sharex=True,
        constrained_layout=True,
        height_ratios=(1, 1),
    )
    top_bars = axes[0].bar(positions, hydropathy, color=bar_colors, width=0.78)
    axes[0].axhline(0, color="#444444", linewidth=0.8)
    axes[0].set(ylabel="Kyte-Doolittle hydropathy")
    axes[0].set_title("Mixed chemistry creates several testable interaction motifs")
    axes[0].spines[["top", "right"]].set_visible(False)
    axes[0].legend(
        handles=[
            Patch(facecolor=color, label=label.replace("_", " "))
            for label, color in chemistry_colors.items()
        ]
        + [
            Patch(
                facecolor="white",
                edgecolor="#111111",
                linewidth=1.8,
                label="top-five candidate",
            )
        ],
        frameon=False,
        fontsize=8,
        ncol=3,
        loc="upper center",
    )
    for bar, position in zip(top_bars, positions):
        if position in candidate_positions:
            bar.set_edgecolor("#111111")
            bar.set_linewidth(1.8)

    bottom_bars = axes[1].bar(positions, confidence, color="#d95f0e", width=0.78)
    axes[1].axhline(70, color="#555555", linestyle="--", linewidth=0.9, label="pLDDT 70")
    axes[1].set(
        xlabel="RAC1B insertion residue",
        ylabel="AlphaFold pLDDT",
        ylim=(0, 100),
        xticks=positions,
        xticklabels=labels,
    )
    axes[1].spines[["top", "right"]].set_visible(False)
    axes[1].legend(frameon=False)
    for label, position in zip(axes[1].get_xticklabels(), positions):
        if position in candidate_positions:
            label.set_fontweight("bold")
    fig.suptitle(
        "RAC1B mutation candidates are chemistry-driven, not fixed-structure claims",
        fontsize=14,
    )
    fig.savefig(figure_dir / "figure4_interface_candidates.png", dpi=200)
    plt.close(fig)

    print(f"Wrote figures to {figure_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
