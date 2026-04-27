#!/usr/bin/env bash
set -euo pipefail

bin=${1:?usage: $0 <phase_mnv_rs_binary>}
repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
fixtures="$repo_root/tests/fixtures"
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

cp "$fixtures/ref.fa" "$tmp/ref.fa"
ref="$tmp/ref.fa"

"$bin" \
  -q \
  -r "$ref" \
  -s S1 \
  --phase-from-bam "$fixtures/read_phase.bam" \
  --no-header \
  "$fixtures/read_phase.vcf" > "$tmp/out.vcf"

diff -u "$fixtures/read_phase.expected.body.vcf" "$tmp/out.vcf"

"$bin" \
  -r "$ref" \
  -s S1 \
  --phase-from-bam "$fixtures/read_phase.bam" \
  --no-header \
  "$fixtures/read_phase.vcf" > "$tmp/out.with-summary.vcf" 2> "$tmp/summary.err"

grep -q 'phase_mnv: bam_phase input=' "$tmp/summary.err"
grep -q 'candidates=4' "$tmp/summary.err"
grep -q 'phased_variants=4' "$tmp/summary.err"
grep -q 'emitted_calls=2' "$tmp/summary.err"

echo "phase_mnv BAM phasing tests passed"
