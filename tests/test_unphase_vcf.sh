#!/usr/bin/env bash
set -euo pipefail

bin=${1:?usage: $0 <unphase_vcf_binary>}
repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
fixtures="$repo_root/tests/fixtures"
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

"$bin" "$fixtures/phased_mnv.vcf" > "$tmp/unphased.vcf"
"$bin" - < "$fixtures/phased_mnv.vcf" > "$tmp/unphased.stdin.vcf"
diff -u "$tmp/unphased.vcf" "$tmp/unphased.stdin.vcf"
! grep -q '##FORMAT=<ID=PS,' "$tmp/unphased.vcf"
! grep -q '##FORMAT=<ID=PQ,' "$tmp/unphased.vcf"
! grep -q '|' "$tmp/unphased.vcf"
grep -q $'\tGT\t0/1$' "$tmp/unphased.vcf"
grep -q $'\tGT\t1/0$' "$tmp/unphased.vcf"
grep -q $'\tGT\t1/1$' "$tmp/unphased.vcf"

"$bin" --keep-phase-tags "$fixtures/phased_mnv.vcf" > "$tmp/keep.vcf"
grep -q '##FORMAT=<ID=PS,' "$tmp/keep.vcf"
grep -q $'\tGT:PS\t0/1:10$' "$tmp/keep.vcf"
! grep -q '|' "$tmp/keep.vcf"

cat > "$tmp/mixed_ploidy.vcf" <<'EOF'
##fileformat=VCFv4.3
##contig=<ID=chr1>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3
chr1	1	.	A	G	.	PASS	.	GT:PS	0|1:10	1:10	.:10
EOF
"$bin" "$tmp/mixed_ploidy.vcf" > "$tmp/mixed_ploidy.out.vcf"
grep -q $'chr1\t1\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\t1\t.' "$tmp/mixed_ploidy.out.vcf"

if command -v bcftools >/dev/null 2>&1; then
  bcftools view -Ob -o "$tmp/phased.bcf" "$fixtures/phased_mnv.vcf"
  "$bin" "$tmp/phased.bcf" > "$tmp/from_bcf.vcf"
  ! grep -q '##FORMAT=<ID=PS,' "$tmp/from_bcf.vcf"
  grep -q $'\tGT\t0/1$' "$tmp/from_bcf.vcf"

  "$bin" -o "$tmp/unphased.vcf.gz" "$fixtures/phased_mnv.vcf"
  bcftools view -H "$tmp/unphased.vcf.gz" > "$tmp/unphased.body"
  grep -q $'\tGT\t0/1$' "$tmp/unphased.body"
fi

echo "unphase_vcf tests passed"
