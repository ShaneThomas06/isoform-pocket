"""Compatibility entry point; implementation lives in isoform_pocket.summarize_p2rank."""

from isoform_pocket.summarize_p2rank import *  # noqa: F403
from isoform_pocket.summarize_p2rank import main

if __name__ == "__main__":
    raise SystemExit(main())
