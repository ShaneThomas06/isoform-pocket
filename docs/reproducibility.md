# Clean-clone verification

Use Python 3.12 (the CI version). From a terminal:

```sh
git clone https://github.com/ShaneThomas06/isoform-pocket.git
cd isoform-pocket
python -m venv .venv
```

Activate with `source .venv/bin/activate` on Linux/macOS or
`.\.venv\Scripts\Activate.ps1` in Windows PowerShell, then run:

```sh
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
python -m pip check
python -m ruff check src scripts tests
python -m pytest -q
python scripts/fetch_inputs.py --help
```

Expected: lint and tests pass, and the final command prints usage and exits
successfully. These checks require no study downloads, GPU, or model weights;
installation requires internet access. Tests use small synthetic fixtures and
temporary outputs. They do not regenerate or overwrite tracked scientific results.
CI runs the same checks on every push and pull request.

This is a software smoke check, not a reproduction of the complete study.
For an exact environment record, save `python -m pip freeze > environment.txt`
alongside your own run outputs, and record the repository commit and input checksums.
Version bounds describe supported dependency families, not an exact lockfile.

## Full analysis

Install `python -m pip install -e ".[pockets,dev]"` for pyKVFinder cavity analysis.
The existing scientific dependency pins are retained. P2Rank additionally needs
the Java runtime installed by `scripts/setup_p2rank.py`; it is not needed for CI.
Follow the README reproduction stages for source downloads and manifest verification.
AlphaFold Server complex inputs must be supplied as documented under
`data/alphafold_server/README.md`; unit tests cannot validate their scientific results.

The `scripts/*.py` commands now delegate to `src/isoform_pocket/`.
Install the project before running them. Reusable helpers can be imported from
`isoform_pocket`, for example `isoform_pocket.compare_structures.kabsch`.
Editable installs retain the checkout-relative data and results paths. With a wheel,
run from the checkout root or set `ISOFORM_POCKET_ROOT` to that checkout before
importing the package. Full analysis commands may regenerate results; use a separate
checkout when comparing against the committed reference outputs.
