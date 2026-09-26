"""Compatibility entry point; implementation lives in isoform_pocket.compare_structures."""

from isoform_pocket.compare_structures import *  # noqa: F403
from isoform_pocket.compare_structures import main

if __name__ == "__main__":
    raise SystemExit(main())
