#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
FIXTURE_DIR="${FIXTURE_DIR:-$ROOT/external/tutorial_bwa_small}"
TMPDIR="${TMPDIR:-$ROOT/.tmp}"
BIN="$ROOT/target/release/bwa-mem2-rs"
OUT_DIR="${OUT_DIR:-$ROOT/.tmp/real_bench}"
REF="$FIXTURE_DIR/ecoli_rel606.fasta.gz"
R1="$FIXTURE_DIR/SRR2584863_1.trim.sub.fastq"
R2="$FIXTURE_DIR/SRR2584863_2.trim.sub.fastq"
PREFIX="$OUT_DIR/ecoli_rel606"
READ_LIMIT="${READ_LIMIT:-0}"

mkdir -p "$TMPDIR" "$OUT_DIR"

if [[ ! -f "$REF" || ! -f "$R1" || ! -f "$R2" ]]; then
  echo "fixture files not found under $FIXTURE_DIR" >&2
  echo "run scripts/download-real-world-fixtures.py first" >&2
  exit 1
fi

CARGO_NET_OFFLINE=true TMPDIR="$TMPDIR" cargo build --release --bin bwa-mem2-rs

if [[ ! -f "${PREFIX}.bwt.2bit.64" ]]; then
  echo "[bench] indexing real-world fixture" >&2
  /usr/bin/time -f "index real %e s" "$BIN" index -p "$PREFIX" "$REF"
fi

BENCH_R1="$R1"
BENCH_R2="$R2"
if [[ "$READ_LIMIT" -gt 0 ]]; then
  BENCH_R1="$OUT_DIR/subset_${READ_LIMIT}_1.fq"
  BENCH_R2="$OUT_DIR/subset_${READ_LIMIT}_2.fq"
  head -n $((READ_LIMIT * 4)) "$R1" >"$BENCH_R1"
  head -n $((READ_LIMIT * 4)) "$R2" >"$BENCH_R2"
  echo "[bench] using paired subset of ${READ_LIMIT} reads" >&2
fi

echo "[bench] paired-end mem -t 1" >&2
/usr/bin/time -f "mem t1 real %e s" "$BIN" mem -o "$OUT_DIR/t1.sam" -t 1 "$PREFIX" "$BENCH_R1" "$BENCH_R2"

echo "[bench] paired-end mem -t 2" >&2
/usr/bin/time -f "mem t2 real %e s" "$BIN" mem -o "$OUT_DIR/t2.sam" -t 2 "$PREFIX" "$BENCH_R1" "$BENCH_R2"

tail -n +3 "$OUT_DIR/t1.sam" >"$OUT_DIR/t1.body.sam"
tail -n +3 "$OUT_DIR/t2.sam" >"$OUT_DIR/t2.body.sam"

cmp "$OUT_DIR/t1.body.sam" "$OUT_DIR/t2.body.sam"
echo "[bench] t1/t2 SAM bodies are identical (ignoring @PG CL header differences)" >&2
