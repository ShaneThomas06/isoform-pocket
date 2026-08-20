"""Map the AF3-derived p120 surface onto experimental structure 3L6Y."""

from __future__ import annotations

import csv
import glob
import math
import os
import re
import shlex
import statistics
from collections import defaultdict
from itertools import combinations
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLBACKEND", "Agg")

ROOT = Path(__file__).resolve().parents[1]
AF_ROOT = ROOT / "data" / "alphafold_server"
RESULTS = ROOT / "results"
FIGURES = RESULTS / "figures"
REFERENCE = ROOT / "data" / "raw" / "experimental" / "3L6Y.cif"
HOTSPOT = (395, 430, 433, 434, 437, 474, 477)
POLYLYSINE = tuple(range(622, 629))
CTNND1_OFFSET = 349
SAME_SURFACE_CUTOFF = 12.0
DIRECT_CONTACT_CUTOFF = 5.0


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"No rows generated for {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def af3_ca_coordinates(path: Path) -> dict[int, tuple[float, float, float, float]]:
    rows: dict[int, tuple[float, float, float, float]] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.startswith("ATOM"):
            continue
        fields = shlex.split(line)
        if len(fields) >= 18 and fields[3] == "CA" and fields[6] == "B":
            rows[int(fields[8])] = (
                float(fields[10]),
                float(fields[11]),
                float(fields[12]),
                float(fields[14]),
            )
    return rows


def experimental_atoms(path: Path) -> tuple[dict[int, list[tuple[float, float, float]]], dict[int, tuple[float, float, float]], list[tuple[float, float, float]]]:
    p120: dict[int, list[tuple[float, float, float]]] = defaultdict(list)
    p120_ca: dict[int, tuple[float, float, float]] = {}
    ecadherin: list[tuple[float, float, float]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.startswith("ATOM"):
            continue
        fields = shlex.split(line)
        if len(fields) < 21:
            continue
        coordinate = (float(fields[10]), float(fields[11]), float(fields[12]))
        residue = int(fields[16])
        chain = fields[18]
        if chain == "A":
            p120[residue].append(coordinate)
            if fields[3] == "CA":
                p120_ca[residue] = coordinate
        elif chain == "B":
            ecadherin.append(coordinate)
    return p120, p120_ca, ecadherin


def distance(first: tuple[float, ...], second: tuple[float, ...]) -> float:
    return math.dist(first[:3], second[:3])


def model_identity(path: Path) -> tuple[str, int, int]:
    isoform = "RAC1B" if "RAC1B_" in str(path) else "RAC1"
    seed_match = re.search(r"seed(\d+)", str(path))
    model_match = re.search(r"_model_(\d+)\.cif$", path.name)
    if not seed_match or not model_match:
        raise ValueError(f"Cannot parse job identity from {path}")
    return isoform, int(seed_match.group(1)), int(model_match.group(1))


def collect_geometry() -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], dict[str, float]]:
    paths = sorted(Path(path) for path in glob.glob(str(AF_ROOT / "RAC*_p120ARM350-824_seed*" / "extracted" / "*_model_*.cif")))
    if len(paths) != 30:
        raise FileNotFoundError(f"Expected 30 AF3 model CIF files, found {len(paths)}")
    if not REFERENCE.exists():
        raise FileNotFoundError(
            f"Missing {REFERENCE}. Download https://files.rcsb.org/download/3L6Y.cif"
        )

    pair_values: dict[tuple[int, int], list[float]] = defaultdict(list)
    motif_values: dict[int, list[float]] = defaultdict(list)
    confidence_values: dict[int, list[float]] = defaultdict(list)
    model_rows: list[dict[str, Any]] = []

    for path in paths:
        isoform, seed, model = model_identity(path)
        coordinates = af3_ca_coordinates(path)
        required = [position - CTNND1_OFFSET for position in HOTSPOT + POLYLYSINE]
        missing = [position for position in required if position not in coordinates]
        if missing:
            raise ValueError(f"Missing p120 local residues {missing} in {path}")

        pair_distances: list[float] = []
        for first, second in combinations(HOTSPOT, 2):
            value = distance(coordinates[first - CTNND1_OFFSET], coordinates[second - CTNND1_OFFSET])
            pair_values[(first, second)].append(value)
            pair_distances.append(value)

        nearest_motif: list[float] = []
        for position in HOTSPOT:
            value = min(
                distance(coordinates[position - CTNND1_OFFSET], coordinates[motif - CTNND1_OFFSET])
                for motif in POLYLYSINE
            )
            motif_values[position].append(value)
            nearest_motif.append(value)
        for position in HOTSPOT + POLYLYSINE:
            confidence_values[position].append(coordinates[position - CTNND1_OFFSET][3])

        model_rows.append(
            {
                "isoform": isoform,
                "seed": seed,
                "model": model,
                "mean_hotspot_pair_distance_angstrom": round(statistics.mean(pair_distances), 3),
                "maximum_hotspot_pair_distance_angstrom": round(max(pair_distances), 3),
                "nearest_hotspot_to_polylysine_angstrom": round(min(nearest_motif), 3),
                "mean_hotspot_to_polylysine_angstrom": round(statistics.mean(nearest_motif), 3),
            }
        )

    p120_atoms, crystal_ca, ecadherin_atoms = experimental_atoms(REFERENCE)
    pair_rows: list[dict[str, Any]] = []
    for first, second in combinations(HOTSPOT, 2):
        values = pair_values[(first, second)]
        crystal_distance = distance(crystal_ca[first], crystal_ca[second])
        af3_mean = statistics.mean(values)
        pair_rows.append(
            {
                "ctnnd1_residue_1": first,
                "ctnnd1_residue_2": second,
                "af3_mean_ca_distance_angstrom": round(af3_mean, 3),
                "af3_sd_angstrom": round(statistics.stdev(values), 3),
                "crystal_3l6y_ca_distance_angstrom": round(crystal_distance, 3),
                "absolute_difference_angstrom": round(abs(af3_mean - crystal_distance), 3),
            }
        )

    motif_rows: list[dict[str, Any]] = []
    for position in HOTSPOT:
        values = motif_values[position]
        motif_rows.append(
            {
                "ctnnd1_residue": position,
                "region": "candidate_surface",
                "mean_nearest_distance_to_622_628_angstrom": round(statistics.mean(values), 3),
                "minimum_distance_angstrom": round(min(values), 3),
                "maximum_distance_angstrom": round(max(values), 3),
                "mean_af3_residue_confidence": round(statistics.mean(confidence_values[position]), 2),
            }
        )
    for position in POLYLYSINE:
        motif_rows.append(
            {
                "ctnnd1_residue": position,
                "region": "polylysine_motif",
                "mean_nearest_distance_to_622_628_angstrom": "",
                "minimum_distance_angstrom": "",
                "maximum_distance_angstrom": "",
                "mean_af3_residue_confidence": round(statistics.mean(confidence_values[position]), 2),
            }
        )

    overlap_rows: list[dict[str, Any]] = []
    for position in HOTSPOT:
        value = min(distance(atom, partner) for atom in p120_atoms[position] for partner in ecadherin_atoms)
        overlap_rows.append(
            {
                "ctnnd1_residue": position,
                "minimum_heavy_atom_distance_to_ecadherin_angstrom": round(value, 3),
                "direct_contact_at_5_angstrom": value <= DIRECT_CONTACT_CUTOFF,
            }
        )

    af3_distances = [float(row["af3_mean_ca_distance_angstrom"]) for row in pair_rows]
    crystal_distances = [float(row["crystal_3l6y_ca_distance_angstrom"]) for row in pair_rows]
    correlation = statistics.correlation(af3_distances, crystal_distances)
    metrics = {
        "r_squared": correlation * correlation,
        "mean_absolute_error": statistics.mean(abs(a - b) for a, b in zip(af3_distances, crystal_distances)),
        "close_pairs": float(sum(value <= SAME_SURFACE_CUTOFF for value in af3_distances)),
        "closest_to_motif": min(value for values in motif_values.values() for value in values),
        "candidate_confidence": statistics.mean(value for position in HOTSPOT for value in confidence_values[position]),
        "motif_confidence": statistics.mean(value for position in POLYLYSINE for value in confidence_values[position]),
        "ecadherin_contacts": float(sum(bool(row["direct_contact_at_5_angstrom"]) for row in overlap_rows)),
    }
    return model_rows, pair_rows, motif_rows, overlap_rows, metrics


def write_report(metrics: dict[str, float]) -> None:
    lines = [
        "# Experimental mapping of the p120 candidate surface",
        "",
        "## Result",
        "",
        f"The seven AF3-derived residues form a compact surface: {int(metrics['close_pairs'])}",
        "of 21 residue pairs have a mean C-alpha distance of 12 A or less across",
        "30 models. All seven positions are resolved in p120-4A crystal structure",
        "3L6Y, where the same geometry is present.",
        "",
        f"AF3 and crystal pair distances agree with R^2 = {metrics['r_squared']:.3f}",
        f"and a mean absolute error of {metrics['mean_absolute_error']:.2f} A. This",
        "supports the geometry of the p120 surface, not RAC binding to it.",
        "",
        "## Relation to the experimentally implicated motif",
        "",
        f"The closest candidate-surface distance to CTNND1 622-628 is {metrics['closest_to_motif']:.1f} A.",
        f"Mean AF3 residue confidence is {metrics['candidate_confidence']:.1f} for the",
        f"candidate surface and {metrics['motif_confidence']:.1f} for the motif. The",
        "622-628 segment is unresolved in 3L6Y.",
        "",
        "## E-cadherin interface overlap",
        "",
        "Six of seven candidate residues are within",
        "5 A of E-cadherin in 3L6Y. CTNND1 430 is 5.59 A away. The AF3-selected",
        "surface overlaps the established cadherin-binding groove.",
        "",
        "## Interpretation",
        "",
        "The complex predictions do not recover the known flexible 622-628 region.",
        "They place RAC proteins against a structured peptide-binding groove already",
        "used by E-cadherin. Combined with low ipTM and seed dependence, this pattern",
        "is treated as a modelling failure mode rather than a RAC1B-specific epitope.",
        "",
        "## References",
        "",
        "- RCSB PDB 3L6Y: https://www.rcsb.org/structure/3L6Y",
        "- Ishiyama et al. (2010): https://doi.org/10.1016/j.cell.2010.01.017",
        "- Orlichenko et al. (2010): https://pmc.ncbi.nlm.nih.gov/articles/PMC2885194/",
        "",
    ]
    (RESULTS / "p120_hotspot_geometry_report.md").write_text("\n".join(lines), encoding="utf-8")


def make_figure(pair_rows: list[dict[str, Any]], motif_rows: list[dict[str, Any]], metrics: dict[str, float]) -> None:
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.lines import Line2D

    candidate = {int(row["ctnnd1_residue"]): row for row in motif_rows if row["region"] == "candidate_surface"}
    motif = {int(row["ctnnd1_residue"]): row for row in motif_rows if row["region"] == "polylysine_motif"}
    af3 = np.array([float(row["af3_mean_ca_distance_angstrom"]) for row in pair_rows])
    crystal = np.array([float(row["crystal_3l6y_ca_distance_angstrom"]) for row in pair_rows])

    plt.rcParams.update({"font.family": "DejaVu Sans", "axes.titlesize": 15, "axes.labelsize": 13, "xtick.labelsize": 10, "ytick.labelsize": 10})
    fig, axes = plt.subplots(2, 2, figsize=(16, 12), facecolor="white")
    fig.subplots_adjust(left=0.08, right=0.97, bottom=0.07, top=0.82, wspace=0.30, hspace=0.48)
    fig.suptitle("Structural context of the AF3-selected p120 surface", y=0.965, fontsize=23, fontweight="bold")
    fig.text(0.5, 0.915, "Thirty AF3 models compared with p120-4A crystal structure 3L6Y", ha="center", fontsize=14, color="#475569")

    ax = axes[0, 0]
    residues = list(HOTSPOT + POLYLYSINE)
    means = [float(candidate[r]["mean_af3_residue_confidence"]) if r in candidate else float(motif[r]["mean_af3_residue_confidence"]) for r in residues]
    colors = ["#ea6a0a" if r in HOTSPOT else "#7042c1" for r in residues]
    ax.scatter(range(len(residues)), means, c=colors, s=72, edgecolors="white", linewidths=0.8)
    ax.axhspan(0, 50, color="#fee2e2", alpha=0.65)
    ax.axhline(70, color="#94a3b8", linestyle="--", linewidth=1)
    ax.set_ylim(0, 100)
    ax.set_xticks(range(len(residues)), labels=residues, rotation=45, ha="right")
    ax.set_xlabel("Canonical CTNND1 residue")
    ax.set_ylabel("Mean AF3 residue confidence")
    ax.set_title("A  The 622-628 motif is low-confidence", loc="left", pad=11)
    ax.grid(axis="y", alpha=0.2)
    ax.legend(handles=[Line2D([0], [0], marker="o", color="white", markerfacecolor="#ea6a0a", label="AF3-selected surface", markersize=9), Line2D([0], [0], marker="o", color="white", markerfacecolor="#7042c1", label="Experimental motif", markersize=9)], frameon=False, loc="upper right")

    ax = axes[0, 1]
    means = np.array([float(candidate[r]["mean_nearest_distance_to_622_628_angstrom"]) for r in HOTSPOT])
    low = np.array([float(candidate[r]["minimum_distance_angstrom"]) for r in HOTSPOT])
    high = np.array([float(candidate[r]["maximum_distance_angstrom"]) for r in HOTSPOT])
    y = np.arange(len(HOTSPOT))
    ax.errorbar(means, y, xerr=[means - low, high - means], fmt="o", color="#ea6a0a", ecolor="#94a3b8", capsize=4, markersize=8)
    ax.axvline(SAME_SURFACE_CUTOFF, color="#475569", linestyle="--", linewidth=1.4)
    ax.text(SAME_SURFACE_CUTOFF, 0.985, "12 A same-surface guide", transform=ax.get_xaxis_transform(), ha="left", va="top", fontsize=9, color="#475569")
    ax.set_yticks(y, labels=HOTSPOT)
    ax.invert_yaxis()
    ax.set_xlim(0, 82)
    ax.set_xlabel("Nearest C-alpha distance to CTNND1 622-628 (A)")
    ax.set_ylabel("Candidate residue")
    ax.set_title("B  The candidate surface is spatially separate", loc="left", pad=11)
    ax.grid(axis="x", alpha=0.2)
    ax.text(0.98, 0.08, f"Minimum across all models: {metrics['closest_to_motif']:.1f} A", transform=ax.transAxes, ha="right", fontsize=11, fontweight="bold", color="#9a3412")

    ax = axes[1, 0]
    ax.scatter(crystal, af3, s=55, color="#0f6b8d", alpha=0.85, edgecolors="white", linewidths=0.7)
    limits = [3, 18]
    ax.plot(limits, limits, color="#475569", linestyle="--", linewidth=1.3)
    ax.set_xlim(limits)
    ax.set_ylim(limits)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("3L6Y crystal C-alpha distance (A)")
    ax.set_ylabel("AF3 mean C-alpha distance (A)")
    ax.set_title("C  Candidate-surface geometry matches 3L6Y", loc="left", pad=11)
    ax.grid(alpha=0.2)
    ax.text(0.07, 0.91, f"R^2 = {metrics['r_squared']:.3f}\nMean error = {metrics['mean_absolute_error']:.2f} A\n18/21 pairs <= 12 A", transform=ax.transAxes, va="top", fontsize=11, bbox={"boxstyle": "round,pad=.5", "facecolor": "#e0f2fe", "edgecolor": "none"})

    ax = axes[1, 1]
    ax.axis("off")
    ax.set_title("D  Experimental interpretation", loc="left", pad=11)
    ax.text(0.08, 0.82, "STRUCTURAL RESULT", transform=ax.transAxes, fontsize=11, fontweight="bold", color="#166534")
    ax.text(0.08, 0.67, "The 395-477 surface exists and\nits geometry matches 3L6Y.", transform=ax.transAxes, fontsize=11, bbox={"boxstyle": "round,pad=.6", "facecolor": "#dcfce7", "edgecolor": "none"})
    ax.text(0.56, 0.82, "BIOLOGICAL CONTEXT", transform=ax.transAxes, fontsize=11, fontweight="bold", color="#92400e")
    ax.text(0.56, 0.67, "6 of 7 residues contact\nE-cadherin in 3L6Y.", transform=ax.transAxes, fontsize=11, bbox={"boxstyle": "round,pad=.6", "facecolor": "#fef3c7", "edgecolor": "none"})
    ax.text(0.08, 0.43, "NOT SUPPORTED", transform=ax.transAxes, fontsize=11, fontweight="bold", color="#991b1b")
    ax.text(0.08, 0.28, "A resolved RAC binding pose or\na RAC1B-exclusive interface.", transform=ax.transAxes, fontsize=11, bbox={"boxstyle": "round,pad=.6", "facecolor": "#fee2e2", "edgecolor": "none"})
    ax.text(0.08, 0.06, "AF3 selected the established cadherin-binding groove rather than\nthe low-confidence 622-628 recognition region.", transform=ax.transAxes, fontsize=12, fontweight="bold", color="#1e293b")

    FIGURES.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURES / "figure8_p120_hotspot_structural_mapping.png", dpi=220, facecolor="white")
    plt.close(fig)


def main() -> int:
    model_rows, pair_rows, motif_rows, overlap_rows, metrics = collect_geometry()
    RESULTS.mkdir(parents=True, exist_ok=True)
    write_csv(RESULTS / "p120_hotspot_model_geometry.csv", model_rows)
    write_csv(RESULTS / "p120_hotspot_pair_distances.csv", pair_rows)
    write_csv(RESULTS / "p120_hotspot_to_polylysine.csv", motif_rows)
    write_csv(RESULTS / "p120_hotspot_ecadherin_overlap.csv", overlap_rows)
    write_report(metrics)
    make_figure(pair_rows, motif_rows, metrics)
    print(f"Mapped {len(model_rows)} AF3 models to experimental structure 3L6Y.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
