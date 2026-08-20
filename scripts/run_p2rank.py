"""Run pinned P2Rank with its AlphaFold-specific prediction profile."""

from __future__ import annotations

import json
import subprocess
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]


def file_map(target: dict[str, Any]) -> dict[str, Path]:
    return {item["kind"]: ROOT / item["path"] for item in target["files"]}


def main() -> int:
    p2rank_home = ROOT / "tools" / "p2rank"
    java_home = ROOT / "tools" / "java17"
    java = java_home / "bin" / "java.exe"
    p2rank_jar = p2rank_home / "bin" / "p2rank.jar"
    if not p2rank_jar.exists() or not java.exists():
        raise RuntimeError("Run scripts/setup_p2rank.py first")

    manifest = json.loads((ROOT / "data" / "manifest.json").read_text(encoding="utf-8"))
    output_root = ROOT / "results" / "p2rank"
    output_root.mkdir(parents=True, exist_ok=True)

    for target in manifest["targets"]:
        label = f'{target["gene"]}_{target["isoform_name"]}'
        output = output_root / label
        output.mkdir(parents=True, exist_ok=True)
        pdb = file_map(target)["alphafold_pdb"]
        command = [
            str(java),
            "-Xmx2048m",
            "-cp",
            f"{p2rank_jar};{p2rank_home / 'bin' / 'lib' / '*'}",
            "cz.siret.prank.program.Main",
            "predict",
            "-f",
            str(pdb),
            "-o",
            str(output),
            "-c",
            "alphafold",
        ]
        print(f"Running P2Rank on {label}", flush=True)
        subprocess.run(command, cwd=p2rank_home, check=True)
    print(f"P2Rank outputs written to {output_root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
