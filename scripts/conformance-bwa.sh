#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
FIXTURE_DIR="${FIXTURE_DIR:-$ROOT/external/tutorial_bwa_small}"
OUT_DIR="${OUT_DIR:-$ROOT/.tmp/conformance_bwa}"
TMPDIR="${TMPDIR:-$ROOT/.tmp}"
RUST_BIN="${RUST_BIN:-$ROOT/target/release/bwa-mem2-rs}"
CPP_BIN="${CPP_BIN:-$ROOT/bwa-mem2/bwa-mem2.avx512bw}"
READ_LIMIT="${READ_LIMIT:-2000}"
BENCH_K="${BENCH_K:-20000000}"
STRICT_KNOWN_GAPS="${STRICT_KNOWN_GAPS:-0}"

REF="$FIXTURE_DIR/ecoli_rel606.fasta.gz"
R1="$FIXTURE_DIR/SRR2584863_1.trim.sub.fastq"
R2="$FIXTURE_DIR/SRR2584863_2.trim.sub.fastq"
PREFIX="$OUT_DIR/ecoli_rel606"
SUB_R1="$OUT_DIR/subset_${READ_LIMIT}_1.fq"
SUB_R2="$OUT_DIR/subset_${READ_LIMIT}_2.fq"
REPORT="$OUT_DIR/report.tsv"

mkdir -p "$TMPDIR" "$OUT_DIR"

if [[ ! -f "$REF" || ! -f "$R1" || ! -f "$R2" ]]; then
  echo "fixture files not found under $FIXTURE_DIR" >&2
  echo "run scripts/download-real-world-fixtures.py first" >&2
  exit 1
fi

if [[ ! -x "$CPP_BIN" ]]; then
  echo "missing executable upstream C++ binary: $CPP_BIN" >&2
  exit 1
fi

CARGO_NET_OFFLINE=true TMPDIR="$TMPDIR" cargo build --release --bin bwa-mem2-rs

if [[ ! -f "${PREFIX}.bwt.2bit.64" ]]; then
  echo "[conf] indexing fixture with Rust binary: $PREFIX" >&2
  "$RUST_BIN" index -p "$PREFIX" "$REF"
fi

if [[ "$READ_LIMIT" -gt 0 ]]; then
  head -n $((READ_LIMIT * 4)) "$R1" >"$SUB_R1"
  head -n $((READ_LIMIT * 4)) "$R2" >"$SUB_R2"
else
  SUB_R1="$R1"
  SUB_R2="$R2"
fi

body_path() {
  local sam="$1"
  local body="$2"
  python3 - "$sam" "$body" <<'PY'
import sys
src, dst = sys.argv[1], sys.argv[2]
with open(src, "r", encoding="utf-8", errors="strict") as inp, open(dst, "w", encoding="utf-8") as out:
    for line in inp:
        if not line.startswith("@"):
            out.write(line)
PY
}

diff_summary() {
  local expected="$1"
  local actual="$2"
  python3 - "$expected" "$actual" <<'PY'
from itertools import zip_longest
import sys

expected, actual = sys.argv[1], sys.argv[2]
diffs = 0
first = None
with open(expected, encoding="utf-8") as fa, open(actual, encoding="utf-8") as fb:
    for i, (a, b) in enumerate(zip_longest(fa, fb), 1):
        if a != b:
            diffs += 1
            if first is None:
                af = (a or "").rstrip("\n").split("\t")[:6]
                bf = (b or "").rstrip("\n").split("\t")[:6]
                first = f"{i}: cpp={af} rust={bf}"
print(f"{diffs}\t{first or ''}")
PY
}

run_case() {
  local name="$1"
  local mode="$2"
  local known_gap="$3"
  shift 3
  local args=("$@")
  local cpp_sam="$OUT_DIR/${name}.cpp.sam"
  local rust_sam="$OUT_DIR/${name}.rust.sam"
  local cpp_body="$OUT_DIR/${name}.cpp.body.sam"
  local rust_body="$OUT_DIR/${name}.rust.body.sam"
  local cpp_log="$OUT_DIR/${name}.cpp.log"
  local rust_log="$OUT_DIR/${name}.rust.log"

  local reads=("$SUB_R1")
  if [[ "$mode" == "pe" ]]; then
    reads=("$SUB_R1" "$SUB_R2")
  fi

  echo "[conf] $name ($mode): ${args[*]}" >&2
  "$CPP_BIN" mem -o "$cpp_sam" "${args[@]}" "$PREFIX" "${reads[@]}" >"$cpp_log" 2>&1
  "$RUST_BIN" mem -o "$rust_sam" "${args[@]}" "$PREFIX" "${reads[@]}" >"$rust_log" 2>&1
  body_path "$cpp_sam" "$cpp_body"
  body_path "$rust_sam" "$rust_body"

  local summary
  summary="$(diff_summary "$cpp_body" "$rust_body")"
  local diffs="${summary%%$'\t'*}"
  local first="${summary#*$'\t'}"
  local status="PASS"
  if [[ "$diffs" != "0" ]]; then
    if [[ "$known_gap" == "1" && "$STRICT_KNOWN_GAPS" != "1" ]]; then
      status="KNOWN-GAP"
    else
      status="FAIL"
    fi
  fi
  printf "%s\t%s\t%s\t%s\t%s\n" "$status" "$name" "$mode" "$diffs" "$first" | tee -a "$REPORT"
}

printf "status\tcase\tmode\tdiff_lines\tfirst_diff\n" >"$REPORT"

failures=0
run_and_count() {
  run_case "$@"
  local status
  status="$(tail -n 1 "$REPORT" | cut -f1)"
  if [[ "$status" == "FAIL" ]]; then
    failures=$((failures + 1))
  fi
}

run_and_count default_pe_t1 pe 0 -t 1 -K "$BENCH_K"
run_and_count default_pe_t2 pe 0 -t 2 -K "$BENCH_K"
run_and_count default_pe_t4 pe 0 -t 4 -K "$BENCH_K"
run_and_count default_se_t1 se 0 -t 1 -K "$BENCH_K"
run_and_count default_se_t2 se 0 -t 2 -K "$BENCH_K"
run_and_count small_batch_pe pe 0 -t 2 -K 50000
run_and_count scoring_A2_pe pe 1 -t 1 -K "$BENCH_K" -A 2
run_and_count scoring_A4_pe pe 1 -t 1 -K "$BENCH_K" -A 4
run_and_count scoring_A5_pe pe 1 -t 1 -K "$BENCH_K" -A 5
run_and_count mismatch_B3_pe pe 1 -t 1 -K "$BENCH_K" -B 3
run_and_count mismatch_B6_pe pe 1 -t 1 -K "$BENCH_K" -B 6
run_and_count gap_open_O4_pe pe 1 -t 1 -K "$BENCH_K" -O 4
run_and_count gap_open_O8_pe pe 1 -t 1 -K "$BENCH_K" -O 8
run_and_count gap_ext_E2_pe pe 1 -t 1 -K "$BENCH_K" -E 2
run_and_count combo_A2_B6_O8_E2_pe pe 1 -t 1 -K "$BENCH_K" -A 2 -B 6 -O 8 -E 2
run_and_count zdrop_150_pe pe 1 -t 1 -K "$BENCH_K" -d 150
run_and_count zdrop_200_pe pe 1 -t 1 -K "$BENCH_K" -d 200

echo "[conf] report: $REPORT" >&2
if [[ "$failures" -ne 0 ]]; then
  echo "[conf] unexpected failures: $failures" >&2
  exit 1
fi
