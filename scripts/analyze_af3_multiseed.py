"""Compatibility entry point; implementation lives in isoform_pocket.analyze_af3_multiseed."""

from isoform_pocket.analyze_af3_multiseed import *  # noqa: F403
from isoform_pocket.analyze_af3_multiseed import main

if __name__ == "__main__":
    raise SystemExit(main())
