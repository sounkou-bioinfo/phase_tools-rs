#!/usr/bin/env bash
set -euo pipefail

bin=${1:?usage: $0 <phase_mnv_rs_binary>}
repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
fixtures="$repo_root/tests/fixtures"
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

if ! command -v bcftools >/dev/null 2>&1; then
  echo "bcftools not found; skipping bcftools norm test"
  exit 0
fi

run_and_check() {
  local name=$1
  shift
  local out="$tmp/$name.vcf"
  local norm="$tmp/$name.norm.vcf"
  local err="$tmp/$name.norm.err"
  local body="$tmp/$name.body"
  local norm_body="$tmp/$name.norm.body"

  "$bin" -q -r "$fixtures/ref.fa" -o "$out" "$@"
  bcftools norm -f "$fixtures/ref.fa" -c x -Ov "$out" > "$norm" 2> "$err"

  grep -v '^#' "$out" > "$body"
  grep -v '^#' "$norm" > "$norm_body"
  diff -u "$body" "$norm_body"

  # The body comparison above is the core assertion. Also require bcftools to
  # report zero realignment and zero mismatch-removal events for this fixture.
  grep -Eq 'realigned/mismatch_removed/.+:\s*[0-9]+/[0-9]+/[0-9]+/0/0/' "$err"
}

run_and_check phased_mnv "$fixtures/phased_mnv.vcf"
run_and_check gap_max1 --max-gap 1 "$fixtures/gap.vcf"
run_and_check complex "$fixtures/complex.vcf"
run_and_check multiallelic "$fixtures/multiallelic.vcf"
run_and_check symbolic_max1 --max-gap 1 "$fixtures/symbolic.vcf"
run_and_check n_base "$fixtures/n_base.vcf"
run_and_check read_phase_bam \
  --sample S1 \
  --phase-from-bam "$fixtures/read_phase.bam" \
  "$fixtures/read_phase.vcf"

if "$bin" --help 2>&1 | grep -q -- '--mnv-algorithm'; then
  run_and_check nirvana_codon \
    --mnv-algorithm nirvana-codon \
    --codon-map "$fixtures/nirvana_codon.codons.tsv" \
    "$fixtures/nirvana_codon.vcf"
fi

echo "phase_mnv bcftools norm tests passed"
