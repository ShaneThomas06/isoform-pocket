"""Resolve study data separately from installed package files."""

import os
from pathlib import Path

# Editable installs preserve the historical checkout-relative behavior.
# Wheel users can set ISOFORM_POCKET_ROOT or run from the study checkout.
_source_root = Path(__file__).resolve().parents[2]
ROOT = Path(os.environ.get("ISOFORM_POCKET_ROOT", str(
    _source_root if (_source_root / "config").is_dir() else Path.cwd()
))).resolve()
