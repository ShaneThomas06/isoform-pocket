"""Compatibility entry point; implementation lives in isoform_pocket.analyze_interface_features."""

from isoform_pocket.analyze_interface_features import *  # noqa: F403
from isoform_pocket.analyze_interface_features import main

if __name__ == "__main__":
    raise SystemExit(main())
