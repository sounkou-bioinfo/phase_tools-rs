#!/usr/bin/env bash
set -euo pipefail

bin=${1:?usage: $0 <phase_mnv_binary>}
repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
fixtures="$repo_root/tests/fixtures"
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

cp "$fixtures/ref.fa" "$tmp/ref.fa"
ref="$tmp/ref.fa"

help="$tmp/help.txt"
"$bin" --help > "$help"
grep -q -- "--unsupported-alleles MODE" "$help"
grep -q -- "--write-index" "$help"
grep -q -- "--no-write-index" "$help"
grep -q -- "combined" "$help"
grep -q -- "--warn-on-n" "$help"
grep -q "Multi-allelic input sites use the ALT allele selected" "$help"
grep -q "unselected ALTs are ignored and output" "$help"
grep -q "Symbolic, breakend, spanning-deletion" "$help"
grep -q "currently not barriers" "$help"
has_nirvana_codon=0
if grep -q -- "--mnv-algorithm MODE" "$help"; then
  has_nirvana_codon=1
  grep -q "nirvana-codon" "$help"
fi
grep -q "output=stdout for VCF stdout" "$help"

run_body() {
  local input=$1
  local output=$2
  shift 2
  "$bin" -q -r "$ref" -s S1 "$@" "$input" > "$tmp/out.vcf"
  grep -v '^#' "$tmp/out.vcf" > "$output" || true
}

run_body "$fixtures/phased_mnv.vcf" "$tmp/phased.body"
diff -u "$fixtures/phased_mnv.expected.body.vcf" "$tmp/phased.body"

run_body "$fixtures/phased_mnv.vcf" "$tmp/combined.body" --emit combined
diff -u "$fixtures/combined.expected.body.vcf" "$tmp/combined.body"

run_body "$fixtures/combined_multiallelic_partial.vcf" "$tmp/combined_multiallelic_partial.body" --emit combined
diff -u "$fixtures/combined_multiallelic_partial.expected.body.vcf" "$tmp/combined_multiallelic_partial.body"

run_body "$fixtures/combined_same_position_records.vcf" "$tmp/combined_same_position_records.body" --emit combined
diff -u "$fixtures/combined_same_position_records.expected.body.vcf" "$tmp/combined_same_position_records.body"

run_body "$fixtures/combined_multisample.vcf" "$tmp/combined_multisample.body" --emit combined
diff -u "$fixtures/combined_multisample.expected.body.vcf" "$tmp/combined_multisample.body"

run_body "$fixtures/combined_multisample_info.vcf" "$tmp/combined_multisample_info.body" --emit combined
diff -u "$fixtures/combined_multisample_info.expected.body.vcf" "$tmp/combined_multisample_info.body"

run_body "$fixtures/combined_residual_other_sample_only.vcf" "$tmp/combined_residual_other_sample_only.body" --emit combined
diff -u "$fixtures/combined_residual_other_sample_only.expected.body.vcf" "$tmp/combined_residual_other_sample_only.body"

run_body "$fixtures/combined_single_sample_info.vcf" "$tmp/combined_single_sample_info.body" --emit combined
diff -u "$fixtures/combined_single_sample_info.expected.body.vcf" "$tmp/combined_single_sample_info.body"

"$bin" -q -r "$ref" -s S1 --emit combined "$fixtures/combined_multisample.vcf" > "$tmp/combined_multisample.vcf"
grep -qx $'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1' "$tmp/combined_multisample.vcf"
! grep -q $'\tS2' "$tmp/combined_multisample.vcf"

"$bin" -q -r "$ref" -s S2 --emit combined "$fixtures/combined_multisample.vcf" > "$tmp/combined_multisample_s2.vcf"
grep -qx $'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS2' "$tmp/combined_multisample_s2.vcf"
grep -v '^#' "$tmp/combined_multisample_s2.vcf" > "$tmp/combined_multisample_s2.body" || true
diff -u "$fixtures/combined_multisample_s2.expected.body.vcf" "$tmp/combined_multisample_s2.body"

if command -v bcftools >/dev/null 2>&1; then
  bcftools view -Ob -o "$tmp/phased.bcf" "$fixtures/phased_mnv.vcf"
  run_body "$tmp/phased.bcf" "$tmp/phased_bcf.body"
  diff -u "$fixtures/phased_mnv.expected.body.vcf" "$tmp/phased_bcf.body"

  bcftools view -Ob -o "$tmp/combined_multiallelic_partial.bcf" "$fixtures/combined_multiallelic_partial.vcf"
  run_body "$tmp/combined_multiallelic_partial.bcf" "$tmp/combined_multiallelic_partial_bcf.body" --emit combined
  diff -u "$fixtures/combined_multiallelic_partial.expected.body.vcf" "$tmp/combined_multiallelic_partial_bcf.body"

  bcftools view -Ob -o "$tmp/combined_multisample_info.bcf" "$fixtures/combined_multisample_info.vcf"
  run_body "$tmp/combined_multisample_info.bcf" "$tmp/combined_multisample_info_bcf.body" --emit combined
  diff -u "$fixtures/combined_multisample_info.expected.body.vcf" "$tmp/combined_multisample_info_bcf.body"
fi

run_body "$fixtures/gap.vcf" "$tmp/gap0.body"
diff -u "$fixtures/gap.max0.expected.body.vcf" "$tmp/gap0.body"

run_body "$fixtures/gap.vcf" "$tmp/gap1.body" --max-gap 1
diff -u "$fixtures/gap.max1.expected.body.vcf" "$tmp/gap1.body"

run_body "$fixtures/complex.vcf" "$tmp/complex.body"
diff -u "$fixtures/complex.expected.body.vcf" "$tmp/complex.body"

run_body "$fixtures/multiallelic.vcf" "$tmp/multiallelic.body"
diff -u "$fixtures/multiallelic.expected.body.vcf" "$tmp/multiallelic.body"

run_body "$fixtures/symbolic.vcf" "$tmp/symbolic.body" --max-gap 1
diff -u "$fixtures/symbolic.max1.expected.body.vcf" "$tmp/symbolic.body"

run_body "$fixtures/n_base.vcf" "$tmp/n_base.body" --warn-on-n 2> "$tmp/n_base.err"
diff -u "$fixtures/n_base.expected.body.vcf" "$tmp/n_base.body"
grep -q "warning: N base in selected allele at chr1:2 hap=2 REF=C ALT=N" "$tmp/n_base.err"

if [[ "$has_nirvana_codon" -eq 1 ]]; then
  run_body "$fixtures/nirvana_codon.vcf" "$tmp/nirvana_codon.body" \
    --mnv-algorithm nirvana-codon \
    --codon-map "$fixtures/nirvana_codon.codons.tsv"
  diff -u "$fixtures/nirvana_codon.expected.body.vcf" "$tmp/nirvana_codon.body"
fi

echo "phase_mnv tests passed"
