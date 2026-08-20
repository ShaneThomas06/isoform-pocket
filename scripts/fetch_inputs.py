"""Download versioned inputs from UniProt, AlphaFold DB, and RCSB PDB."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
USER_AGENT = "IsoformPocket/0.1 (reproducible academic project)"


def fetch_bytes(url: str, timeout: int = 90) -> tuple[bytes, dict[str, str]]:
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        content = response.read()
        headers = {key.lower(): value for key, value in response.headers.items()}
    return content, headers


def checksum(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def save_download(
    *, content: bytes, path: Path, url: str, kind: str, headers: dict[str, str]
) -> dict[str, Any]:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(content)
    return {
        "kind": kind,
        "path": path.relative_to(ROOT).as_posix(),
        "url": url,
        "sha256": checksum(content),
        "bytes": len(content),
        "last_modified": headers.get("last-modified"),
    }


def select_alphafold_record(records: list[dict[str, Any]], query: str) -> dict[str, Any]:
    if not records:
        raise ValueError(f"AlphaFold DB returned no record for {query}")
    exact = [
        record
        for record in records
        if record.get("uniprotAccession") == query
        or record.get("entryId", "").startswith(f"AF-{query}-")
    ]
    return exact[0] if exact else records[0]


def fetch_target(target: dict[str, str], raw_dir: Path) -> dict[str, Any]:
    accession = target["uniprot_accession"]
    query = target["alphafold_query"]
    destination = raw_dir / target["gene"] / target["isoform_name"]
    files: list[dict[str, Any]] = []

    fasta_url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
    content, headers = fetch_bytes(fasta_url)
    files.append(
        save_download(
            content=content,
            path=destination / f"{accession}.fasta",
            url=fasta_url,
            kind="uniprot_fasta",
            headers=headers,
        )
    )

    api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{query}"
    api_content, api_headers = fetch_bytes(api_url)
    record = select_alphafold_record(json.loads(api_content), query)
    normalized = (json.dumps(record, indent=2, sort_keys=True) + "\n").encode()
    files.append(
        save_download(
            content=normalized,
            path=destination / f"{query}.alphafold.json",
            url=api_url,
            kind="alphafold_metadata",
            headers=api_headers,
        )
    )

    urls = {
        "alphafold_pdb": record.get("pdbUrl"),
        "alphafold_cif": record.get("cifUrl"),
        "alphafold_pae": record.get("paeDocUrl"),
    }
    suffixes = {
        "alphafold_pdb": ".pdb",
        "alphafold_cif": ".cif",
        "alphafold_pae": ".pae.json",
    }
    for kind, url in urls.items():
        if not url:
            raise ValueError(f"{record.get('entryId', query)} has no {kind} URL")
        content, headers = fetch_bytes(url)
        files.append(
            save_download(
                content=content,
                path=destination / f"{query}{suffixes[kind]}",
                url=url,
                kind=kind,
                headers=headers,
            )
        )

    return {**target, "alphafold_entry_id": record.get("entryId"), "files": files}


def fetch_experimental(row: dict[str, str], raw_dir: Path) -> dict[str, Any]:
    pdb_id = row["pdb_id"].upper()
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    content, headers = fetch_bytes(url)
    path = raw_dir / "experimental" / f"{pdb_id}.pdb"
    file_record = save_download(
        content=content, path=path, url=url, kind="experimental_pdb", headers=headers
    )
    return {**row, "pdb_id": pdb_id, "files": [file_record]}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-dir", type=Path, default=ROOT / "data" / "raw")
    parser.add_argument("--manifest", type=Path, default=ROOT / "data" / "manifest.json")
    args = parser.parse_args()

    manifest: dict[str, Any] = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "targets": [],
        "experimental_structures": [],
        "errors": [],
    }

    for target in read_csv(ROOT / "config" / "targets.csv"):
        label = f'{target["gene"]}/{target["isoform_name"]}'
        print(f"Fetching {label}", flush=True)
        try:
            manifest["targets"].append(fetch_target(target, args.raw_dir))
        except (OSError, ValueError, json.JSONDecodeError, urllib.error.URLError) as exc:
            manifest["errors"].append({"item": label, "error": str(exc)})
            print(f"ERROR {label}: {exc}", file=sys.stderr)

    for row in read_csv(ROOT / "config" / "experimental_structures.csv"):
        label = f'PDB/{row["pdb_id"]}'
        print(f"Fetching {label}", flush=True)
        try:
            manifest["experimental_structures"].append(fetch_experimental(row, args.raw_dir))
        except (OSError, ValueError, urllib.error.URLError) as exc:
            manifest["errors"].append({"item": label, "error": str(exc)})
            print(f"ERROR {label}: {exc}", file=sys.stderr)

    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {args.manifest}")
    return 1 if manifest["errors"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
