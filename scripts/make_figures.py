"""Compatibility entry point; implementation lives in isoform_pocket.make_figures."""

from isoform_pocket.make_figures import *  # noqa: F403
from isoform_pocket.make_figures import main

if __name__ == "__main__":
    raise SystemExit(main())
