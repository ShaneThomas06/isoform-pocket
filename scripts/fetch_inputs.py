"""Compatibility entry point; implementation lives in isoform_pocket.fetch_inputs."""

from isoform_pocket.fetch_inputs import *  # noqa: F403
from isoform_pocket.fetch_inputs import main

if __name__ == "__main__":
    raise SystemExit(main())
