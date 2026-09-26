"""Compatibility entry point; implementation lives in isoform_pocket.verify_manifest."""

from isoform_pocket.verify_manifest import *  # noqa: F403
from isoform_pocket.verify_manifest import main

if __name__ == "__main__":
    raise SystemExit(main())
