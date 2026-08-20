"""Download pinned, portable P2Rank and Temurin Java distributions."""

from __future__ import annotations

import hashlib
import shutil
import tarfile
import urllib.request
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools"
DOWNLOADS = TOOLS / ".downloads"

P2RANK_URL = "https://github.com/rdk/p2rank/releases/download/2.5.1/p2rank_2.5.1.tar.gz"
P2RANK_SHA256 = "d243f2d9036ac053fefb9407b5fe1c85f4fe077c519fd975ac585e995feab274"
JRE_URL = (
    "https://github.com/adoptium/temurin17-binaries/releases/download/"
    "jdk-17.0.19%2B10/OpenJDK17U-jre_x64_windows_hotspot_17.0.19_10.zip"
)
JRE_SHA256 = "79a598e1fbb4e16582d92c4ee22280a3c4d72fd52606e1e46b1223c0fe53b0da"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def download(url: str, destination: Path, expected: str) -> None:
    if destination.exists() and sha256(destination) == expected:
        print(f"Using verified {destination.name}")
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(url, headers={"User-Agent": "IsoformPocket/0.1"})
    print(f"Downloading {destination.name}")
    with urllib.request.urlopen(request, timeout=180) as response, destination.open("wb") as output:
        shutil.copyfileobj(response, output, length=1024 * 1024)
    observed = sha256(destination)
    if observed != expected:
        destination.unlink(missing_ok=True)
        raise RuntimeError(f"Checksum mismatch for {destination.name}: {observed}")


def main() -> int:
    p2rank_archive = DOWNLOADS / "p2rank_2.5.1.tar.gz"
    jre_archive = DOWNLOADS / "OpenJDK17U-jre_x64_windows_hotspot_17.0.19_10.zip"
    download(P2RANK_URL, p2rank_archive, P2RANK_SHA256)
    download(JRE_URL, jre_archive, JRE_SHA256)

    p2rank_destination = TOOLS / "p2rank"
    java_destination = TOOLS / "java17"
    if not p2rank_destination.exists():
        print("Extracting P2Rank")
        with tarfile.open(p2rank_archive, "r:gz") as archive:
            archive.extractall(TOOLS, filter="data")
        extracted = TOOLS / "p2rank_2.5.1"
        extracted.rename(p2rank_destination)
    if not java_destination.exists():
        print("Extracting Java")
        with zipfile.ZipFile(jre_archive) as archive:
            archive.extractall(TOOLS)
        extracted = next(TOOLS.glob("jdk-17.0.19+10-jre*"))
        extracted.rename(java_destination)

    java = java_destination / "bin" / "java.exe"
    if not java.exists() or not (p2rank_destination / "bin" / "p2rank.jar").exists():
        raise RuntimeError("Portable tool extraction did not create expected files")
    print(f"P2Rank ready at {p2rank_destination}")
    print(f"Java ready at {java}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
