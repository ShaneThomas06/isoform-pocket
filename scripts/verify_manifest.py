"""Verify downloaded inputs against the recorded SHA-256 manifest."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "data" / "manifest.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def iter_records(manifest: dict):
    for target in manifest.get("targets", []):
        yield from target.get("files", [])
    for structure in manifest.get("experimental_structures", []):
        yield from structure.get("files", [])


def main() -> int:
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    failures: list[str] = []
    records = list(iter_records(manifest))

    for record in records:
        relative = Path(record["path"])
        path = ROOT / relative
        if not path.is_file():
            failures.append(f"MISSING {relative}")
            continue
        actual_size = path.stat().st_size
        if actual_size != record["bytes"]:
            failures.append(
                f"SIZE {relative}: expected {record['bytes']}, got {actual_size}"
            )
            continue
        actual_hash = sha256(path)
        if actual_hash != record["sha256"]:
            failures.append(
                f"SHA256 {relative}: expected {record['sha256']}, got {actual_hash}"
            )

    if failures:
        print("\n".join(failures))
        print(f"Verification failed: {len(failures)} of {len(records)} files.")
        return 1

    print(f"Verified {len(records)} files against data/manifest.json.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
