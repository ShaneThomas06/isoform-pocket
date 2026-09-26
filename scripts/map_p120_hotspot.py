"""Compatibility entry point; implementation lives in isoform_pocket.map_p120_hotspot."""

from isoform_pocket.map_p120_hotspot import *  # noqa: F403
from isoform_pocket.map_p120_hotspot import main

if __name__ == "__main__":
    raise SystemExit(main())
