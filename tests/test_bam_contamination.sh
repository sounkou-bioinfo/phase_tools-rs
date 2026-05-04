#!/usr/bin/env bash
set -euo pipefail

bin=${1:?usage: $0 <bam_contamination_binary>}
repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
fixtures="$repo_root/tests/fixtures"
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

cat > "$tmp/anchors.tsv" <<'EOF'
chrom	pos	ref	alt	gt	ref_af
chr1	1	A	G	1/1	0.5
chr1	2	C	T	1|1	0.5
chr1	3	G	T	0/1	0.5
EOF

"$bin" \
  --reference "$fixtures/ref.fa" \
  --bam "$fixtures/read_phase.bam" \
  --anchors "$tmp/anchors.tsv" > "$tmp/out.tsv"

grep -qx $'#anchors\t3' "$tmp/out.tsv"
grep -qx $'chrom\tpos\tref\talt\tgt\tclass\tref_af\tobservations\tref_count\talt_count\tother_count\tignored_count\tref_fraction\talt_fraction\tother_fraction\tcharr_like_component' "$tmp/out.tsv"
grep -qx $'chr1\t1\tA\tG\t1/1\thom_alt\t0.500000\t4\t2\t2\t0\t0\t0.500000\t0.500000\t0.000000\t1.000000' "$tmp/out.tsv"
grep -qx $'chr1\t2\tC\tT\t1|1\thom_alt\t0.500000\t4\t2\t2\t0\t0\t0.500000\t0.500000\t0.000000\t1.000000' "$tmp/out.tsv"
grep -qx $'chr1\t3\tG\tT\t0/1\thet\t0.500000\t4\t4\t0\t0\t0\t1.000000\t0.000000\t0.000000\tNA' "$tmp/out.tsv"
grep -qx $'#hom_alt_sites\t2' "$tmp/out.tsv"
grep -qx $'#hom_alt_observations\t8' "$tmp/out.tsv"
grep -qx $'#hom_alt_mean_ref_balance\t0.500000' "$tmp/out.tsv"
grep -qx $'#charr_like_sites\t2' "$tmp/out.tsv"
grep -qx $'#charr_like_mean\t1.000000' "$tmp/out.tsv"

cat > "$tmp/anchors.vcf" <<'EOF'
##fileformat=VCFv4.3
##contig=<ID=chr1,length=12>
##INFO=<ID=REF_AF,Number=1,Type=Float,Description="Reference allele frequency">
##INFO=<ID=AF,Number=A,Type=Float,Description="Alternate allele frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
chr1	1	.	A	G	.	.	REF_AF=0.5	GT	1/1
chr1	2	.	C	T	.	.	AF=0.5	GT	1|1
chr1	3	.	G	T	.	.	AF=0.5	GT	0/1
chr1	4	.	T	A	.	.	REF_AF=.	GT	1/1
EOF
"$bin" \
  --reference "$fixtures/ref.fa" \
  --bam "$fixtures/read_phase.bam" \
  --anchors "$tmp/anchors.vcf" \
  --sample SAMPLE > "$tmp/vcf.out"
grep -qx $'#anchors\t4' "$tmp/vcf.out"
grep -qx $'chr1\t1\tA\tG\t1/1\thom_alt\t0.500000\t4\t2\t2\t0\t0\t0.500000\t0.500000\t0.000000\t1.000000' "$tmp/vcf.out"
grep -qx $'chr1\t2\tC\tT\t1|1\thom_alt\t0.500000\t4\t2\t2\t0\t0\t0.500000\t0.500000\t0.000000\t1.000000' "$tmp/vcf.out"
grep -qx $'chr1\t3\tG\tT\t0/1\thet\t0.500000\t4\t4\t0\t0\t0\t1.000000\t0.000000\t0.000000\tNA' "$tmp/vcf.out"
grep -qx $'chr1\t4\tT\tA\t1/1\thom_alt\tNA\t4\t2\t0\t2\t0\t1.000000\t0.000000\t0.500000\tNA' "$tmp/vcf.out"

gzip -c "$tmp/anchors.vcf" > "$tmp/anchors.vcf.bgz"
"$bin" \
  --reference "$fixtures/ref.fa" \
  --bam "$fixtures/read_phase.bam" \
  --anchors "$tmp/anchors.vcf.bgz" \
  --sample SAMPLE > "$tmp/vcf_bgz.out"
grep -qx $'#anchors\t4' "$tmp/vcf_bgz.out"
grep -qx $'chr1\t4\tT\tA\t1/1\thom_alt\tNA\t4\t2\t0\t2\t0\t1.000000\t0.000000\t0.500000\tNA' "$tmp/vcf_bgz.out"

cat > "$tmp/ref_mismatch.vcf" <<'EOF'
##fileformat=VCFv4.3
##contig=<ID=chr1,length=12>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
chr1	1	.	C	G	.	.	.	GT	1/1
EOF
if "$bin" --reference "$fixtures/ref.fa" --bam "$fixtures/read_phase.bam" --anchors "$tmp/ref_mismatch.vcf" > "$tmp/ref_mismatch_vcf.out" 2> "$tmp/ref_mismatch_vcf.err"; then
  echo "bam_contamination unexpectedly accepted VCF anchor/FASTA REF mismatch" >&2
  exit 1
fi
grep -q 'anchor REF mismatch' "$tmp/ref_mismatch_vcf.err"

"$bin" \
  --reference "$fixtures/ref.fa" \
  --bam "$fixtures/read_phase.bam" \
  --anchors "$tmp/anchors.tsv" \
  --min-baseq 41 > "$tmp/baseq.tsv"
grep -qx $'chr1\t1\tA\tG\t1/1\thom_alt\t0.500000\t0\t0\t0\t0\t4\tNA\tNA\tNA\tNA' "$tmp/baseq.tsv"
grep -qx $'#hom_alt_observations\t0' "$tmp/baseq.tsv"
grep -qx $'#charr_like_mean\tNA' "$tmp/baseq.tsv"

"$bin" \
  --reference "$fixtures/ref.fa" \
  --bam "$fixtures/read_phase.bam" \
  --anchors "$tmp/anchors.tsv" \
  --min-mapq 61 > "$tmp/mapq.tsv"
grep -qx $'chr1\t1\tA\tG\t1/1\thom_alt\t0.500000\t0\t0\t0\t0\t0\tNA\tNA\tNA\tNA' "$tmp/mapq.tsv"

echo -e 'chrom\tpos\tref\talt\tgt\nchr1\t1\tC\tG\t1/1' > "$tmp/ref_mismatch.tsv"
if "$bin" --reference "$fixtures/ref.fa" --bam "$fixtures/read_phase.bam" --anchors "$tmp/ref_mismatch.tsv" > "$tmp/ref_mismatch.out" 2> "$tmp/ref_mismatch.err"; then
  echo "bam_contamination unexpectedly accepted anchor/FASTA REF mismatch" >&2
  exit 1
fi
grep -q 'anchor REF mismatch' "$tmp/ref_mismatch.err"

cat > "$tmp/no_af.tsv" <<'EOF'
chrom	pos	ref	alt	gt
chr1	1	A	G	1/1
EOF
"$bin" --reference "$fixtures/ref.fa" --bam "$fixtures/read_phase.bam" --anchors "$tmp/no_af.tsv" > "$tmp/no_af.out"
grep -qx $'chr1\t1\tA\tG\t1/1\thom_alt\tNA\t4\t2\t2\t0\t0\t0.500000\t0.500000\t0.000000\tNA' "$tmp/no_af.out"

cat > "$tmp/bad.tsv" <<'EOF'
chrom	pos	ref	alt	gt	ref_af
chr1	1	A	G	1/2	0.5
EOF
if "$bin" --reference "$fixtures/ref.fa" --bam "$fixtures/read_phase.bam" --anchors "$tmp/bad.tsv" > "$tmp/bad.out" 2> "$tmp/bad.err"; then
  echo "bam_contamination unexpectedly accepted unsupported genotype" >&2
  exit 1
fi
grep -q 'currently supports only 0 and 1 alleles' "$tmp/bad.err"

echo "bam_contamination tests passed"
