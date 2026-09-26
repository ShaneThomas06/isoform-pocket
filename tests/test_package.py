"""Regression checks for installed-package and legacy script behavior."""

import os
from pathlib import Path
import subprocess
import sys

import isoform_pocket.paths


ROOT = Path(__file__).resolve().parents[1]


def test_editable_install_keeps_checkout_paths():
    assert isoform_pocket.paths.ROOT == ROOT


def test_explicit_data_root_from_another_directory(tmp_path):
    env = {**os.environ, "ISOFORM_POCKET_ROOT": str(tmp_path)}
    result = subprocess.run(
        [sys.executable, "-c", "from isoform_pocket.paths import ROOT; print(ROOT)"],
        cwd=tmp_path, env=env, capture_output=True, text=True, check=True,
    )
    assert Path(result.stdout.strip()) == tmp_path


def test_legacy_script_runs_outside_checkout(tmp_path):
    result = subprocess.run(
        [sys.executable, str(ROOT / "scripts/fetch_inputs.py"), "--help"],
        cwd=tmp_path, capture_output=True, text=True, check=True,
    )
    assert "usage:" in result.stdout.lower()


def test_legacy_helper_import_is_same_function():
    # scripts is deliberately not installed as a generic top-level package.
    sys.path.insert(0, str(ROOT))
    try:
        from scripts.compare_structures import kabsch as legacy
        from isoform_pocket.compare_structures import kabsch
        assert legacy is kabsch
    finally:
        sys.path.pop(0)
