"""Compatibility entry point; implementation lives in isoform_pocket.validate_cavities."""

from isoform_pocket.validate_cavities import *  # noqa: F403
from isoform_pocket.validate_cavities import main

if __name__ == "__main__":
    raise SystemExit(main())
