#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
FIXTURE_DIR="${FIXTURE_DIR:-$ROOT/external/ecoli_miseq_2x250}"
OUT_DIR="${OUT_DIR:-$ROOT/.tmp/regression_simd16_a5}"
TMPDIR="${TMPDIR:-$ROOT/.tmp}"
RUST_BIN="${RUST_BIN:-$ROOT/target/release/bwa-mem2-rs}"
CPP_BIN="${CPP_BIN:-$ROOT/bwa-mem2/bwa-mem2.avx512bw}"
READ_LIMIT="${READ_LIMIT:-2000}"
OPTIONAL_READ_LIMIT="${OPTIONAL_READ_LIMIT:-5000}"
BENCH_K="${BENCH_K:-20000000}"

REF="${REF:-$FIXTURE_DIR/ecoli_rel606.fasta.gz}"
R1="${R1:-$FIXTURE_DIR/SRR13321180_1.sub${OPTIONAL_READ_LIMIT}.fastq}"
R2="${R2:-$FIXTURE_DIR/SRR13321180_2.sub${OPTIONAL_READ_LIMIT}.fastq}"
PREFIX="${PREFIX:-$OUT_DIR/ecoli_miseq_2x250}"
SUB_R1="$OUT_DIR/subset_${READ_LIMIT}_1.fq"
SUB_R2="$OUT_DIR/subset_${READ_LIMIT}_2.fq"

mkdir -p "$TMPDIR" "$OUT_DIR"

if [[ ! -f "$REF" || ! -f "$R1" || ! -f "$R2" ]]; then
  echo "fixture files not found under $FIXTURE_DIR" >&2
  echo "run scripts/download-real-world-fixtures.py --extra miseq-2x250 --extra-read-limit $OPTIONAL_READ_LIMIT first" >&2
  exit 1
fi

if [[ ! -x "$CPP_BIN" ]]; then
  echo "missing executable upstream C++ binary: $CPP_BIN" >&2
  exit 1
fi

CARGO_NET_OFFLINE=true TMPDIR="$TMPDIR" cargo build --release --bin bwa-mem2-rs

if [[ ! -f "${PREFIX}.bwt.2bit.64" ]]; then
  echo "[regression] indexing fixture with Rust binary: $PREFIX" >&2
  "$RUST_BIN" index -p "$PREFIX" "$REF"
fi

head -n $((READ_LIMIT * 4)) "$R1" >"$SUB_R1"
head -n $((READ_LIMIT * 4)) "$R2" >"$SUB_R2"

CPP_SAM="$OUT_DIR/a5.cpp.sam"
RUST_SAM="$OUT_DIR/a5.rust.sam"
CPP_BODY="$OUT_DIR/a5.cpp.body.sam"
RUST_BODY="$OUT_DIR/a5.rust.body.sam"

echo "[regression] MiSeq 2x250 PE -A 5 with default SIMD16 dispatch" >&2
"$CPP_BIN" mem -o "$CPP_SAM" -t 1 -K "$BENCH_K" -A 5 "$PREFIX" "$SUB_R1" "$SUB_R2" >"$OUT_DIR/a5.cpp.log" 2>&1
"$RUST_BIN" mem -o "$RUST_SAM" -t 1 -K "$BENCH_K" -A 5 "$PREFIX" "$SUB_R1" "$SUB_R2" >"$OUT_DIR/a5.rust.log" 2>&1

python3 - "$CPP_SAM" "$CPP_BODY" "$RUST_SAM" "$RUST_BODY" <<'PY'
import sys

for src, dst in ((sys.argv[1], sys.argv[2]), (sys.argv[3], sys.argv[4])):
    with open(src, encoding="utf-8") as inp, open(dst, "w", encoding="utf-8") as out:
        for line in inp:
            if not line.startswith("@"):
                out.write(line)
PY

if ! cmp -s "$CPP_BODY" "$RUST_BODY"; then
  python3 - "$CPP_BODY" "$RUST_BODY" <<'PY'
from itertools import zip_longest
import sys

diffs = 0
first = None
with open(sys.argv[1], encoding="utf-8") as expected, open(sys.argv[2], encoding="utf-8") as actual:
    for i, (a, b) in enumerate(zip_longest(expected, actual), 1):
        if a != b:
            diffs += 1
            if first is None:
                first = (i, (a or "").rstrip("\n").split("\t")[:12], (b or "").rstrip("\n").split("\t")[:12])
print(f"[regression] FAIL diff_lines={diffs}", file=sys.stderr)
print(f"[regression] first_diff={first}", file=sys.stderr)
PY
  exit 1
fi

echo "[regression] PASS MiSeq -A 5 SAM bodies are byte-identical" >&2
