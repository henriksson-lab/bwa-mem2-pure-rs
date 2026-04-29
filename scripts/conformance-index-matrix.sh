#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMPDIR="${TMPDIR:-$ROOT/.tmp}"
OUT_DIR="${OUT_DIR:-$ROOT/.tmp/conformance_index_matrix}"
RUST_BIN="${RUST_BIN:-$ROOT/target/release/bwa-mem2-rs}"
CPP_BIN="${CPP_BIN:-$ROOT/bwa-mem2/bwa-mem2.avx512bw}"
REPORT="$OUT_DIR/report.tsv"

mkdir -p "$TMPDIR" "$OUT_DIR"

if [[ ! -x "$CPP_BIN" ]]; then
  echo "missing executable upstream C++ binary: $CPP_BIN" >&2
  exit 1
fi

CARGO_NET_OFFLINE=true TMPDIR="$TMPDIR" cargo build --release --bin bwa-mem2-rs

write_fixtures() {
  local fixture_dir="$1"
  mkdir -p "$fixture_dir"

  printf ">chr1 simple reference\nACGTACGTACGTACGT\n" >"$fixture_dir/single.fa"

  {
    printf ">chr1 first contig\n"
    printf "ACGTACGTACGT\n"
    printf ">chr2 second contig\n"
    printf "TTTTCCCCAAAAGGGG\n"
    printf ">plasmid third contig\n"
    printf "GATTACAGATTACA\n"
  } >"$fixture_dir/multi.fa"

  {
    printf ">ambiguous with runs\n"
    printf "ACGTNNNNACGTMRWSYKBDHVNACGT\n"
    printf ">lowercase_and_masked\n"
    printf "aaaaccccNNNNggggtttt\n"
  } >"$fixture_dir/ambiguous.fa"

  {
    printf ">one_base\n"
    printf "A\n"
    printf ">two_bases\n"
    printf "CG\n"
    printf ">three_bases\n"
    printf "TTA\n"
    printf ">four_bases\n"
    printf "ACGT\n"
  } >"$fixture_dir/short.fa"
}

compare_case() {
  local name="$1"
  local fasta="$2"
  local rust_prefix="$OUT_DIR/$name/rust/ref"
  local cpp_prefix="$OUT_DIR/$name/cpp/ref"
  local status="PASS"
  local details=""
  local files=(pac bwt.2bit.64 ann amb 0123)

  mkdir -p "$OUT_DIR/$name/rust" "$OUT_DIR/$name/cpp"
  echo "[index-conf] $name" >&2
  "$RUST_BIN" index -p "$rust_prefix" "$fasta" >"$OUT_DIR/$name/rust.index.log" 2>&1
  "$CPP_BIN" index -p "$cpp_prefix" "$fasta" >"$OUT_DIR/$name/cpp.index.log" 2>&1

  for ext in "${files[@]}"; do
    local rust_file="$rust_prefix.$ext"
    local cpp_file="$cpp_prefix.$ext"
    if [[ ! -f "$rust_file" || ! -f "$cpp_file" ]]; then
      status="FAIL"
      details="${details}${ext}:missing;"
      continue
    fi
    if ! cmp -s "$cpp_file" "$rust_file"; then
      status="FAIL"
      local cpp_sum rust_sum
      cpp_sum="$(md5sum "$cpp_file" | cut -d' ' -f1)"
      rust_sum="$(md5sum "$rust_file" | cut -d' ' -f1)"
      details="${details}${ext}:cpp=${cpp_sum},rust=${rust_sum};"
    fi
  done

  printf "%s\t%s\t%s\n" "$status" "$name" "${details:-byte-identical}" | tee -a "$REPORT"
  [[ "$status" == "PASS" ]]
}

FIXTURE_DIR="$OUT_DIR/fixtures"
write_fixtures "$FIXTURE_DIR"

printf "status\tcase\tdetails\n" >"$REPORT"
failures=0

for name in single multi ambiguous short; do
  if ! compare_case "$name" "$FIXTURE_DIR/$name.fa"; then
    failures=$((failures + 1))
  fi
done

echo "[index-conf] report: $REPORT" >&2
if [[ "$failures" -ne 0 ]]; then
  echo "[index-conf] $failures case(s) failed" >&2
  exit 1
fi
