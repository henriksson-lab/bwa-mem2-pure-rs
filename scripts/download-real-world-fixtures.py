#!/usr/bin/env python3
import hashlib
import json
import os
import sys
import tarfile
import urllib.request
from pathlib import Path


REFERENCE_URL = (
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/017/985/"
    "GCA_000017985.1_ASM1798v1/GCA_000017985.1_ASM1798v1_genomic.fna.gz"
)
SUBSET_TAR_URL = "https://ndownloader.figshare.com/files/14418248"
PAIR_FILES = [
    "SRR2584863_1.trim.sub.fastq",
    "SRR2584863_2.trim.sub.fastq",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download(url: str, dest: Path) -> None:
    with urllib.request.urlopen(url) as response, dest.open("wb") as out:
        while True:
            chunk = response.read(1024 * 1024)
            if not chunk:
                break
            out.write(chunk)


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    out_dir = root / "external" / "tutorial_bwa_small"
    out_dir.mkdir(parents=True, exist_ok=True)

    manifest = {}

    ref_path = out_dir / "ecoli_rel606.fasta.gz"
    if not ref_path.exists():
        print(f"downloading {ref_path.name} from {REFERENCE_URL}", file=sys.stderr)
        download(REFERENCE_URL, ref_path)
    manifest[ref_path.name] = {
        "url": REFERENCE_URL,
        "size": ref_path.stat().st_size,
        "sha256": sha256(ref_path),
    }

    tar_path = out_dir / "trimmed_fastq_small.tar.gz"
    if not tar_path.exists():
        print(f"downloading {tar_path.name} from {SUBSET_TAR_URL}", file=sys.stderr)
        download(SUBSET_TAR_URL, tar_path)

    with tarfile.open(tar_path, "r:gz") as tf:
        members = {Path(member.name).name: member for member in tf.getmembers()}
        for pair_name in PAIR_FILES:
            pair_path = out_dir / pair_name
            if not pair_path.exists():
                member = members.get(pair_name)
                if member is None:
                    raise FileNotFoundError(f"{pair_name} not found in {tar_path}")
                print(f"extracting {pair_name} from {tar_path.name}", file=sys.stderr)
                with tf.extractfile(member) as src, pair_path.open("wb") as dst:
                    if src is None:
                        raise FileNotFoundError(f"unable to extract {pair_name} from {tar_path}")
                    while True:
                        chunk = src.read(1024 * 1024)
                        if not chunk:
                            break
                        dst.write(chunk)
            manifest[pair_name] = {
                "url": SUBSET_TAR_URL,
                "size": pair_path.stat().st_size,
                "sha256": sha256(pair_path),
            }
    manifest[tar_path.name] = {
        "url": SUBSET_TAR_URL,
        "size": tar_path.stat().st_size,
        "sha256": sha256(tar_path),
    }

    manifest_path = out_dir / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
