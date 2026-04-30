#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT_DIR="${OUT_DIR:-$ROOT/.tmp/median_bench}"
TMPDIR="${TMPDIR:-$ROOT/.tmp}"
RUST_BIN="${RUST_BIN:-$ROOT/target/release/bwa-mem2-rs}"
CPP_BIN="${CPP_BIN:-$ROOT/bwa-mem2/bwa-mem2.avx512bw}"
RUNS="${RUNS:-3}"
RUN_CPP="${RUN_CPP:-1}"
BENCH_K="${BENCH_K:-20000000}"
OPTIONAL_READ_LIMIT="${OPTIONAL_READ_LIMIT:-5000}"
REPORT="$OUT_DIR/report.tsv"

mkdir -p "$TMPDIR" "$OUT_DIR"

if [[ "$RUNS" -lt 1 ]]; then
  echo "RUNS must be >= 1" >&2
  exit 1
fi
if [[ "$RUN_CPP" == "1" && ! -x "$CPP_BIN" ]]; then
  echo "missing executable upstream C++ binary: $CPP_BIN" >&2
  exit 1
fi

CARGO_NET_OFFLINE=true TMPDIR="$TMPDIR" cargo build --release --bin bwa-mem2-rs

body_path() {
  local sam="$1"
  local body="$2"
  python3 - "$sam" "$body" <<'PY'
import sys
src, dst = sys.argv[1:]
with open(src, "r", encoding="utf-8", errors="strict") as inp, open(dst, "w", encoding="utf-8") as out:
    for line in inp:
        if not line.startswith("@"):
            out.write(line)
PY
}

time_command() {
  local label="$1"
  shift
  local start end
  local log_path="${BENCH_LOG:-/dev/null}"
  start="$(date +%s.%N)"
  "$@" >"$log_path" 2>&1
  end="$(date +%s.%N)"
  python3 - "$label" "$start" "$end" <<'PY'
import sys
label, start, end = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
print(f"{label}\t{end - start:.3f}", flush=True)
PY
}

median() {
  python3 - "$@" <<'PY'
import statistics
import sys
values = [float(x) for x in sys.argv[1:] if x]
print(f"{statistics.median(values):.3f}")
PY
}

ensure_index() {
  local prefix="$1"
  local ref="$2"
  if [[ ! -f "${prefix}.bwt.2bit.64" ]]; then
    echo "[bench] indexing $ref -> $prefix" >&2
    "$RUST_BIN" index -p "$prefix" "$ref"
  fi
}

compare_bodies() {
  local expected="$1"
  local actual="$2"
  python3 - "$expected" "$actual" <<'PY'
from itertools import zip_longest
import sys

expected, actual = sys.argv[1:]
diffs = 0
first = None
with open(expected, encoding="utf-8") as a, open(actual, encoding="utf-8") as b:
    for i, (x, y) in enumerate(zip_longest(a, b), 1):
        if x != y:
            diffs += 1
            if first is None:
                first = i
print(f"{diffs}\t{first or ''}")
PY
}

run_case() {
  local dataset="$1"
  local case_name="$2"
  local prefix="$3"
  local r1="$4"
  local r2="$5"
  shift 5
  local args=("$@")
  local case_dir="$OUT_DIR/$dataset/$case_name"
  mkdir -p "$case_dir"

  local rust_times=()
  local cpp_times=()
  local rust_sam cpp_sam rust_body cpp_body

  echo "[bench] $dataset/$case_name: ${args[*]}" >&2
  for run in $(seq 1 "$RUNS"); do
    rust_sam="$case_dir/rust_${run}.sam"
    rust_times+=("$(BENCH_LOG="$case_dir/rust_${run}.log" time_command "rust" "$RUST_BIN" mem -o "$rust_sam" "${args[@]}" "$prefix" "$r1" "$r2" | cut -f2)")
    echo "[bench] $dataset/$case_name rust run $run/${RUNS}: ${rust_times[-1]}s" >&2
  done

  local diff_lines="NA"
  local first_diff=""
  if [[ "$RUN_CPP" == "1" ]]; then
    for run in $(seq 1 "$RUNS"); do
      cpp_sam="$case_dir/cpp_${run}.sam"
      cpp_times+=("$(BENCH_LOG="$case_dir/cpp_${run}.log" time_command "cpp" "$CPP_BIN" mem -o "$cpp_sam" "${args[@]}" "$prefix" "$r1" "$r2" | cut -f2)")
      echo "[bench] $dataset/$case_name cpp run $run/${RUNS}: ${cpp_times[-1]}s" >&2
    done
    rust_body="$case_dir/rust_1.body.sam"
    cpp_body="$case_dir/cpp_1.body.sam"
    body_path "$case_dir/rust_1.sam" "$rust_body"
    body_path "$case_dir/cpp_1.sam" "$cpp_body"
    local summary
    summary="$(compare_bodies "$cpp_body" "$rust_body")"
    diff_lines="${summary%%$'\t'*}"
    first_diff="${summary#*$'\t'}"
    if [[ "$diff_lines" != "0" ]]; then
      echo "[bench] parity failure for $dataset/$case_name: diff_lines=$diff_lines first=$first_diff" >&2
      return 1
    fi
  fi

  local rust_median cpp_median ratio
  rust_median="$(median "${rust_times[@]}")"
  if [[ "$RUN_CPP" == "1" ]]; then
    cpp_median="$(median "${cpp_times[@]}")"
    ratio="$(python3 - "$rust_median" "$cpp_median" <<'PY'
import sys
rust, cpp = map(float, sys.argv[1:])
print(f"{rust / cpp:.3f}")
PY
)"
  else
    cpp_median="NA"
    ratio="NA"
  fi

  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$dataset" "$case_name" "${args[*]}" "$RUNS" "$rust_median" "$cpp_median" "$ratio" "$diff_lines" "$first_diff" \
    | tee -a "$REPORT"
}

printf "dataset\tcase\targs\truns\trust_median_s\tcpp_median_s\trust_cpp_ratio\tdiff_lines\tfirst_diff\n" >"$REPORT"

tutorial_ref="$ROOT/external/tutorial_bwa_small/ecoli_rel606.fasta.gz"
tutorial_r1="$ROOT/external/tutorial_bwa_small/SRR2584863_1.trim.sub.fastq"
tutorial_r2="$ROOT/external/tutorial_bwa_small/SRR2584863_2.trim.sub.fastq"
tutorial_prefix="$OUT_DIR/index/tutorial_bwa_small"
if [[ -f "$tutorial_ref" && -f "$tutorial_r1" && -f "$tutorial_r2" ]]; then
  mkdir -p "$(dirname "$tutorial_prefix")"
  ensure_index "$tutorial_prefix" "$tutorial_ref"
  run_case tutorial_bwa_small default_pe_t1 "$tutorial_prefix" "$tutorial_r1" "$tutorial_r2" -t 1 -K "$BENCH_K"
else
  echo "[bench] skipping tutorial_bwa_small: missing fixture files" >&2
fi

miseq_ref="$ROOT/external/ecoli_miseq_2x250/ecoli_rel606.fasta.gz"
miseq_r1="$ROOT/external/ecoli_miseq_2x250/SRR13321180_1.sub${OPTIONAL_READ_LIMIT}.fastq"
miseq_r2="$ROOT/external/ecoli_miseq_2x250/SRR13321180_2.sub${OPTIONAL_READ_LIMIT}.fastq"
miseq_prefix="$OUT_DIR/index/ecoli_miseq_2x250"
if [[ -f "$miseq_ref" && -f "$miseq_r1" && -f "$miseq_r2" ]]; then
  mkdir -p "$(dirname "$miseq_prefix")"
  ensure_index "$miseq_prefix" "$miseq_ref"
  run_case ecoli_miseq_2x250 default_pe_t1 "$miseq_prefix" "$miseq_r1" "$miseq_r2" -t 1 -K "$BENCH_K"
  run_case ecoli_miseq_2x250 seed_len_k5_pe "$miseq_prefix" "$miseq_r1" "$miseq_r2" -t 1 -K "$BENCH_K" -k 5
else
  echo "[bench] skipping ecoli_miseq_2x250: missing fixture files" >&2
fi

echo "[bench] report: $REPORT" >&2
