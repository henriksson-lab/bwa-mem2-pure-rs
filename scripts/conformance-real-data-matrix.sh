#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT_ROOT="${OUT_ROOT:-$ROOT/.tmp/conformance_real_data_matrix}"
READ_LIMIT="${READ_LIMIT:-2000}"
OPTIONAL_READ_LIMIT="${OPTIONAL_READ_LIMIT:-5000}"
STRICT_KNOWN_GAPS="${STRICT_KNOWN_GAPS:-1}"

run_dataset() {
  local name="$1"
  local fixture_dir="$2"
  local ref="$3"
  local r1="$4"
  local r2="$5"
  local run_pe="$6"
  local run_se="$7"

  if [[ ! -f "$ref" || ! -f "$r1" ]]; then
    echo "[matrix] skip $name: missing $ref or $r1" >&2
    return 0
  fi
  if [[ "$run_pe" == "1" && ! -f "$r2" ]]; then
    echo "[matrix] skip $name: missing paired read $r2" >&2
    return 0
  fi

  echo "[matrix] dataset $name" >&2
  FIXTURE_DIR="$fixture_dir" \
  REF="$ref" \
  R1="$r1" \
  R2="$r2" \
  PREFIX_NAME="$name" \
  OUT_DIR="$OUT_ROOT/$name" \
  READ_LIMIT="$READ_LIMIT" \
  STRICT_KNOWN_GAPS="$STRICT_KNOWN_GAPS" \
  RUN_PE="$run_pe" \
  RUN_SE="$run_se" \
    "$ROOT/scripts/conformance-bwa.sh"
}

run_dataset \
  tutorial_bwa_small \
  "$ROOT/external/tutorial_bwa_small" \
  "$ROOT/external/tutorial_bwa_small/ecoli_rel606.fasta.gz" \
  "$ROOT/external/tutorial_bwa_small/SRR2584863_1.trim.sub.fastq" \
  "$ROOT/external/tutorial_bwa_small/SRR2584863_2.trim.sub.fastq" \
  1 \
  1

run_dataset \
  ecoli_miseq_2x250 \
  "$ROOT/external/ecoli_miseq_2x250" \
  "$ROOT/external/ecoli_miseq_2x250/ecoli_rel606.fasta.gz" \
  "$ROOT/external/ecoli_miseq_2x250/SRR13321180_1.sub${OPTIONAL_READ_LIMIT}.fastq" \
  "$ROOT/external/ecoli_miseq_2x250/SRR13321180_2.sub${OPTIONAL_READ_LIMIT}.fastq" \
  1 \
  0

echo "[matrix] reports under $OUT_ROOT" >&2
