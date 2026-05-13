#!/usr/bin/env bash
set -euo pipefail

bin=${1:?usage: $0 <phase_mnv_rs_binary>}
repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
fixtures="$repo_root/tests/fixtures"
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

if ! command -v bcftools >/dev/null 2>&1; then
  echo "bcftools not found; skipping compressed VCF/BCF output checks"
  exit 0
fi

cp "$fixtures/ref.fa" "$tmp/ref.fa"
ref="$tmp/ref.fa"

"$bin" -q -r "$ref" -s S1 --threads 2 -o "$tmp/out.vcf.gz" "$fixtures/phased_mnv.vcf"
bcftools index -f "$tmp/out.vcf.gz"
bcftools view -H "$tmp/out.vcf.gz" > "$tmp/out.vcf.gz.body"
diff -u "$fixtures/phased_mnv.expected.body.vcf" "$tmp/out.vcf.gz.body"

"$bin" -q -r "$ref" -s S1 --threads 2 --write-index -o "$tmp/self_index.vcf.gz" "$fixtures/phased_mnv.vcf"
test -s "$tmp/self_index.vcf.gz.csi"
bcftools index --stats "$tmp/self_index.vcf.gz" > "$tmp/self_index.stats"
grep -q $'^chr1\t\.\t3$' "$tmp/self_index.stats"

"$bin" -q -r "$ref" -s S1 --threads 2 --write-index=tbi -o "$tmp/self_index_tbi.vcf.gz" "$fixtures/phased_mnv.vcf"
test -s "$tmp/self_index_tbi.vcf.gz.tbi"
bcftools index --stats "$tmp/self_index_tbi.vcf.gz" > "$tmp/self_index_tbi.stats"
grep -q $'^chr1\t\.\t3$' "$tmp/self_index_tbi.stats"

"$bin" -q -r "$ref" -s S1 --threads 2 -o "$tmp/out.vcf.bgz" "$fixtures/phased_mnv.vcf"
bcftools index -f "$tmp/out.vcf.bgz"
bcftools view -H "$tmp/out.vcf.bgz" > "$tmp/out.vcf.bgz.body"
diff -u "$fixtures/phased_mnv.expected.body.vcf" "$tmp/out.vcf.bgz.body"

"$bin" \
  -q \
  -r "$ref" \
  -s S1 \
  --threads 2 \
  --emit all-sites \
  --phase-from-bam "$fixtures/read_phase.bam" \
  -o "$tmp/all_sites.vcf.gz" \
  "$fixtures/read_phase.vcf"
bcftools index -f "$tmp/all_sites.vcf.gz"
bcftools view -H "$tmp/all_sites.vcf.gz" > "$tmp/all_sites.body"
grep -qx $'chr1\t1\t.\tA\tG\t.\tPASS\t.\tGT:PS\t0|1:1' "$tmp/all_sites.body"
grep -qx $'chr1\t5\t.\tA\tG\t.\tPASS\t.\tGT:PS\t1|0:1' "$tmp/all_sites.body"

"$bin" -q -r "$ref" -s S1 --threads 2 -o "$tmp/out.bcf" "$fixtures/phased_mnv.vcf"
bcftools index -f "$tmp/out.bcf"
bcftools view -H "$tmp/out.bcf" > "$tmp/out.bcf.body"
diff -u "$fixtures/phased_mnv.expected.body.vcf" "$tmp/out.bcf.body"

"$bin" -q -r "$ref" -s S1 --threads 2 --write-index -o "$tmp/self_index.bcf" "$fixtures/phased_mnv.vcf"
test -s "$tmp/self_index.bcf.csi"
bcftools index --stats "$tmp/self_index.bcf" > "$tmp/self_index_bcf.stats"
grep -q $'^chr1\t\.\t3$' "$tmp/self_index_bcf.stats"

"$bin" -q -r "$ref" -s S1 --threads 2 --emit combined --write-index -o "$tmp/combined.vcf.gz" "$fixtures/phased_mnv.vcf"
test -s "$tmp/combined.vcf.gz.csi"
bcftools view -H "$tmp/combined.vcf.gz" > "$tmp/combined.body"
diff -u "$fixtures/combined.expected.body.vcf" "$tmp/combined.body"

"$bin" -r "$ref" -s S1 --threads 2 -o "$tmp/summary.vcf.gz" "$fixtures/phased_mnv.vcf" \
  > "$tmp/summary.stdout" 2> "$tmp/summary.stderr"
test ! -s "$tmp/summary.stdout"
grep -q 'output_format=vcf.gz threads=2' "$tmp/summary.stderr"

echo "phase_mnv output format tests passed"
