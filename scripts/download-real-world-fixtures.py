#!/usr/bin/env python3
import argparse
import gzip
import hashlib
import json
import shutil
import sys
import tarfile
import urllib.parse
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
MISEQ_2X250_RUN = "SRR13321180"


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


def ena_read_run(run_accession: str) -> dict:
    fields = ",".join(
        [
            "run_accession",
            "fastq_ftp",
            "fastq_md5",
            "fastq_bytes",
            "read_count",
            "base_count",
            "library_layout",
            "instrument_platform",
            "instrument_model",
        ]
    )
    query = urllib.parse.urlencode(
        {
            "accession": run_accession,
            "result": "read_run",
            "fields": fields,
            "format": "json",
        }
    )
    url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?{query}"
    with urllib.request.urlopen(url) as response:
        data = json.loads(response.read().decode("utf-8"))
    if not data:
        raise RuntimeError(f"ENA returned no metadata for {run_accession}")
    return data[0]


def stream_fastq_subset_gz(url: str, dest: Path, read_limit: int) -> None:
    records = 0
    with urllib.request.urlopen(url) as response, gzip.GzipFile(fileobj=response) as gz, dest.open(
        "wb"
    ) as out:
        while records < read_limit:
            record = [gz.readline() for _ in range(4)]
            if not record[0]:
                break
            if any(line == b"" for line in record):
                raise EOFError(f"truncated FASTQ record while reading {url}")
            for line in record:
                out.write(line)
            records += 1
    if records == 0:
        raise RuntimeError(f"downloaded no FASTQ records from {url}")


def download_tutorial_fixture(root: Path) -> Path:
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
    return out_dir


def download_miseq_2x250_fixture(root: Path, read_limit: int) -> Path:
    tutorial_dir = download_tutorial_fixture(root)
    out_dir = root / "external" / "ecoli_miseq_2x250"
    out_dir.mkdir(parents=True, exist_ok=True)

    ref_src = tutorial_dir / "ecoli_rel606.fasta.gz"
    ref_path = out_dir / "ecoli_rel606.fasta.gz"
    if not ref_path.exists():
        shutil.copyfile(ref_src, ref_path)

    metadata = ena_read_run(MISEQ_2X250_RUN)
    fastq_urls = [f"https://{url}" for url in metadata["fastq_ftp"].split(";") if url]
    if len(fastq_urls) != 2:
        raise RuntimeError(f"expected two FASTQ URLs for {MISEQ_2X250_RUN}, got {fastq_urls}")

    r1 = out_dir / f"{MISEQ_2X250_RUN}_1.sub{read_limit}.fastq"
    r2 = out_dir / f"{MISEQ_2X250_RUN}_2.sub{read_limit}.fastq"
    for url, dest in zip(fastq_urls, [r1, r2]):
        if not dest.exists():
            print(f"streaming first {read_limit} reads from {url}", file=sys.stderr)
            stream_fastq_subset_gz(url, dest, read_limit)

    manifest = {
        ref_path.name: {
            "source": str(ref_src),
            "size": ref_path.stat().st_size,
            "sha256": sha256(ref_path),
        },
        r1.name: {
            "url": fastq_urls[0],
            "run_accession": MISEQ_2X250_RUN,
            "records": read_limit,
            "size": r1.stat().st_size,
            "sha256": sha256(r1),
        },
        r2.name: {
            "url": fastq_urls[1],
            "run_accession": MISEQ_2X250_RUN,
            "records": read_limit,
            "size": r2.stat().st_size,
            "sha256": sha256(r2),
        },
        "ena_metadata": metadata,
    }
    manifest_path = out_dir / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return out_dir


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--extra",
        choices=["miseq-2x250"],
        action="append",
        default=[],
        help="download an additional optional real-data fixture",
    )
    parser.add_argument(
        "--extra-read-limit",
        type=int,
        default=5000,
        help="number of read pairs to keep for optional FASTQ fixtures",
    )
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[1]
    out_dir = download_tutorial_fixture(root)
    extra_dirs = []
    if "miseq-2x250" in args.extra:
        extra_dirs.append(download_miseq_2x250_fixture(root, args.extra_read_limit))

    print(out_dir)
    for path in extra_dirs:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
