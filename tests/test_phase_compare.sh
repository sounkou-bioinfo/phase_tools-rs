#!/usr/bin/env bash
set -euo pipefail

bin=${1:?usage: $0 <phase_compare_binary>}
repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
fixtures="$repo_root/tests/fixtures"
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

"$bin" \
  --sample S1 \
  -r "$fixtures/ref.fa" \
  --engine xcmp \
  --no-roc \
  --no-decompose \
  --switch-error-bed "$tmp/switches.bed" \
  --tsv-pairwise "$tmp/pairs.tsv" \
  --summary-tsv "$tmp/summary.copy.tsv" \
  "$fixtures/phase_compare_truth.vcf" \
  "$fixtures/phase_compare_query.vcf" > "$tmp/summary.tsv"

cmp "$tmp/summary.tsv" "$tmp/summary.copy.tsv"

"$bin" \
  --sample S1 \
  -o "$tmp/prefix" \
  "$fixtures/phase_compare_truth.vcf" \
  "$fixtures/phase_compare_query.vcf" > "$tmp/prefix.stdout.tsv"
test -s "$tmp/prefix.summary.tsv"

header=$(head -1 "$tmp/summary.tsv")
[[ "$header" == chrom$'\t'truth_records* ]]

total=$(awk -F'\t' '$1=="TOTAL" {print}' "$tmp/summary.tsv")
[[ -n "$total" ]]

# Columns:
# 1 chrom, 2 truth_records, 3 query_records, 4 common_records,
# 7 common_het, 11 both_phased_het_with_ps, 12 intersection_blocks,
# 13 intersection_variants, 14 assessed_pairs, 15 phase_match_pairs,
# 16 switch_errors, 17 switch_rate, 18 blockwise_hamming, 19 hamming_rate.
# Singleton PS blocks are counted as phased sites but not as intersection blocks.
awk -F'\t' '$1=="TOTAL" {
  if ($2 != 5 || $3 != 5 || $4 != 5 || $7 != 5 || $11 != 4 ||
      $12 != 1 || $13 != 3 || $14 != 2 || $15 != 1 || $16 != 1 ||
      $17 != "0.500000" || $18 != 1 || $19 != "0.333333") exit 1
}' "$tmp/summary.tsv"

wc -l < "$tmp/switches.bed" | grep -qx '1'
grep -qx $'chr1\t1\t3' "$tmp/switches.bed"

wc -l < "$tmp/pairs.tsv" | grep -qx '3'
grep -q $'chr1\t2\t3\t1\t1\t0\t1\tswitch' "$tmp/pairs.tsv"

echo "phase_compare tests passed"
