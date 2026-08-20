# Raw input data

Large source files are excluded from version control. The analysis scripts expect the following local inputs.

## Downloaded by scripts

Run `python scripts/fetch_inputs.py` from the repository root. Source URLs and SHA-256 checksums are recorded in `data/manifest.json`.

## Experimental p120 structure

Download the official RCSB mmCIF entry for 3L6Y:

```powershell
New-Item -ItemType Directory -Force data/raw/experimental
Invoke-WebRequest -Uri "https://files.rcsb.org/download/3L6Y.cif" -OutFile "data/raw/experimental/3L6Y.cif"
```

Expected path: `data/raw/experimental/3L6Y.cif`.

Downloaded structures retain the licensing and usage terms of their source databases.
