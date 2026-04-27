#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
FIXTURE_DIR="${FIXTURE_DIR:-$ROOT/external/tutorial_bwa_small}"
TMPDIR="${TMPDIR:-$ROOT/.tmp}"
BIN="$ROOT/target/release/bwa-mem2-rs"
CPP_BIN="${CPP_BIN:-$ROOT/bwa-mem2/bwa-mem2.avx512bw}"
OUT_DIR="${OUT_DIR:-$ROOT/.tmp/real_bench}"
REF="$FIXTURE_DIR/ecoli_rel606.fasta.gz"
R1="$FIXTURE_DIR/SRR2584863_1.trim.sub.fastq"
R2="$FIXTURE_DIR/SRR2584863_2.trim.sub.fastq"
PREFIX="$OUT_DIR/ecoli_rel606"
READ_LIMIT="${READ_LIMIT:-0}"
BENCH_K="${BENCH_K:-20000000}"
RUN_CPP="${RUN_CPP:-1}"

mkdir -p "$TMPDIR" "$OUT_DIR"

if [[ ! -f "$REF" || ! -f "$R1" || ! -f "$R2" ]]; then
  echo "fixture files not found under $FIXTURE_DIR" >&2
  echo "run scripts/download-real-world-fixtures.py first" >&2
  exit 1
fi

if [[ "$RUN_CPP" == "1" && ! -x "$CPP_BIN" ]]; then
  echo "missing executable upstream C++ binary: $CPP_BIN" >&2
  exit 1
fi

run_timed() {
  local label="$1"
  shift
  local start end
  start="$(date +%s.%N)"
  "$@"
  end="$(date +%s.%N)"
  python3 - "$label" "$start" "$end" <<'PY'
import sys
label, start, end = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
elapsed = end - start
print(f"{label}\t{elapsed:.2f}", flush=True)
PY
}

body_path() {
  local sam="$1"
  local body="$2"
  python3 - "$sam" "$body" <<'PY'
import sys
src, dst = sys.argv[1], sys.argv[2]
with open(src, encoding="utf-8") as inp, open(dst, "w", encoding="utf-8") as out:
    for line in inp:
        if not line.startswith("@"):
            out.write(line)
PY
}

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

echo "[bench] paired-end mem -t 1 (-K $BENCH_K)" >&2
rust_t1="$(run_timed "rust_t1" "$BIN" mem -o "$OUT_DIR/t1.sam" -t 1 -K "$BENCH_K" "$PREFIX" "$BENCH_R1" "$BENCH_R2" | tee "$OUT_DIR/rust_t1.time.tsv" | cut -f2)"

echo "[bench] paired-end mem -t 2 (-K $BENCH_K)" >&2
rust_t2="$(run_timed "rust_t2" "$BIN" mem -o "$OUT_DIR/t2.sam" -t 2 -K "$BENCH_K" "$PREFIX" "$BENCH_R1" "$BENCH_R2" | tee "$OUT_DIR/rust_t2.time.tsv" | cut -f2)"

body_path "$OUT_DIR/t1.sam" "$OUT_DIR/t1.body.sam"
body_path "$OUT_DIR/t2.sam" "$OUT_DIR/t2.body.sam"

cmp "$OUT_DIR/t1.body.sam" "$OUT_DIR/t2.body.sam"
echo "[bench] t1/t2 SAM bodies are identical (ignoring @PG CL header differences)" >&2

if [[ "$RUN_CPP" == "1" ]]; then
  echo "[bench] upstream C++ paired-end mem -t 1 (-K $BENCH_K)" >&2
  cpp_t1="$(run_timed "cpp_t1" "$CPP_BIN" mem -o "$OUT_DIR/cpp_t1.sam" -t 1 -K "$BENCH_K" "$PREFIX" "$BENCH_R1" "$BENCH_R2" | tee "$OUT_DIR/cpp_t1.time.tsv" | cut -f2)"

  echo "[bench] upstream C++ paired-end mem -t 2 (-K $BENCH_K)" >&2
  cpp_t2="$(run_timed "cpp_t2" "$CPP_BIN" mem -o "$OUT_DIR/cpp_t2.sam" -t 2 -K "$BENCH_K" "$PREFIX" "$BENCH_R1" "$BENCH_R2" | tee "$OUT_DIR/cpp_t2.time.tsv" | cut -f2)"

  body_path "$OUT_DIR/cpp_t1.sam" "$OUT_DIR/cpp_t1.body.sam"
  body_path "$OUT_DIR/cpp_t2.sam" "$OUT_DIR/cpp_t2.body.sam"
  cmp "$OUT_DIR/t1.body.sam" "$OUT_DIR/cpp_t1.body.sam"
  cmp "$OUT_DIR/t2.body.sam" "$OUT_DIR/cpp_t2.body.sam"
  echo "[bench] Rust and upstream C++ SAM bodies are byte-identical" >&2

  python3 - "$rust_t1" "$cpp_t1" "$rust_t2" "$cpp_t2" "$OUT_DIR/report.tsv" <<'PY'
import sys
rust_t1, cpp_t1, rust_t2, cpp_t2 = map(float, sys.argv[1:5])
report = sys.argv[5]
rows = [
    ("mem -t 1", rust_t1, cpp_t1),
    ("mem -t 2", rust_t2, cpp_t2),
]
with open(report, "w", encoding="utf-8") as out:
    out.write("command\trust_seconds\tcpp_seconds\trust_vs_cpp\n")
    for command, rust, cpp in rows:
        out.write(f"{command}\t{rust:.2f}\t{cpp:.2f}\t{rust / cpp:.2f}x\n")
for command, rust, cpp in rows:
    print(f"[bench] {command}: Rust {rust:.2f}s, C++ {cpp:.2f}s, Rust/C++ {rust / cpp:.2f}x", file=sys.stderr)
PY
  echo "[bench] report: $OUT_DIR/report.tsv" >&2
fi
