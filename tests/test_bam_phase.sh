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
  -q \
  -r "$ref" \
  -s S1 \
  --phase-from-bam "$fixtures/read_phase.bam" \
  --phase-algorithm greedy \
  --no-header \
  "$fixtures/read_phase.vcf" > "$tmp/out.greedy.vcf"
diff -u "$fixtures/read_phase.expected.body.vcf" "$tmp/out.greedy.vcf"

if command -v samtools >/dev/null 2>&1; then

cat > "$tmp/blank_bridge.vcf" <<'VCF'
##fileformat=VCFv4.3
##contig=<ID=chr1>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
chr1	1	.	A	C	.	PASS	.	GT	0/1
chr1	2	.	C	T	.	PASS	.	GT	0/1
chr1	3	.	G	T	.	PASS	.	GT	0/1
VCF
cat > "$tmp/blank_bridge.sam" <<'SAM'
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:12
@RG	ID:S1	SM:S1
blank_h1_a	0	chr1	1	60	1M1D1M	*	0	0	AT	II	RG:Z:S1
blank_h2_a	0	chr1	1	60	1M1D1M	*	0	0	CG	II	RG:Z:S1
blank_h1_b	0	chr1	2	60	2M	*	0	0	TT	II	RG:Z:S1
blank_h2_b	0	chr1	2	60	2M	*	0	0	CG	II	RG:Z:S1
SAM
samtools view -b "$tmp/blank_bridge.sam" | samtools sort -o "$tmp/blank_bridge.bam"
samtools index "$tmp/blank_bridge.bam"
"$bin" \
  -q \
  -r "$ref" \
  -s S1 \
  --phase-from-bam "$tmp/blank_bridge.bam" \
  --phase-internal-downsampling 23 \
  --emit all-sites \
  "$tmp/blank_bridge.vcf" > "$tmp/blank_bridge.out.vcf"
grep -v '^#' "$tmp/blank_bridge.out.vcf" > "$tmp/blank_bridge.body.vcf"
cat > "$tmp/blank_bridge.expected.vcf" <<'VCF'
chr1	1	.	A	C	.	PASS	.	GT:PS	0|1:1
chr1	2	.	C	T	.	PASS	.	GT:PS	1|0:1
chr1	3	.	G	T	.	PASS	.	GT:PS	1|0:1
VCF
diff -u "$tmp/blank_bridge.expected.vcf" "$tmp/blank_bridge.body.vcf"

"$bin" \
  -q \
  -r "$ref" \
  -s S1 \
  --phase-from-bam "$tmp/blank_bridge.bam" \
  --tag HP \
  --only-snvs \
  --mapq 0 \
  --ignore-read-groups \
  --output-read-list "$tmp/selected_reads.tsv" \
  --emit all-sites \
  "$tmp/blank_bridge.vcf" > "$tmp/blank_bridge.hp.vcf"
grep -v '^#' "$tmp/blank_bridge.hp.vcf" > "$tmp/blank_bridge.hp.body.vcf"
cat > "$tmp/blank_bridge.hp.expected.vcf" <<'VCF'
chr1	1	.	A	C	.	PASS	.	GT:HP	0|1:1-1,1-2
chr1	2	.	C	T	.	PASS	.	GT:HP	1|0:1-2,1-1
chr1	3	.	G	T	.	PASS	.	GT:HP	1|0:1-2,1-1
VCF
diff -u "$tmp/blank_bridge.hp.expected.vcf" "$tmp/blank_bridge.hp.body.vcf"
test -s "$tmp/selected_reads.tsv"
grep -q '^#read_name' "$tmp/selected_reads.tsv"

python3 - <<'PY' > "$tmp/poly.fa"
print('>chrR')
print('A' * 30)
PY
samtools faidx "$tmp/poly.fa"
cat > "$tmp/leftshifted_insertion.vcf" <<'VCF'
##fileformat=VCFv4.3
##contig=<ID=chrR>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
chrR	2	.	A	C	.	PASS	.	GT	0/1
chrR	20	.	A	AA	.	PASS	.	GT	0/1
VCF
python3 - <<'PY' > "$tmp/leftshifted_insertion.sam"
print('@HD\tVN:1.6\tSO:coordinate')
print('@SQ\tSN:chrR\tLN:30')
print('@RG\tID:S1\tSM:S1')
print('hap_alt_leftshifted\t0\tchrR\t1\t60\t19M1I11M\t*\t0\t0\t' + 'AC' + 'A' * 29 + '\t' + 'I' * 31 + '\tRG:Z:S1')
print('hap_ref\t0\tchrR\t1\t60\t30M\t*\t0\t0\t' + 'A' * 30 + '\t' + 'I' * 30 + '\tRG:Z:S1')
PY
samtools view -b "$tmp/leftshifted_insertion.sam" | samtools sort -o "$tmp/leftshifted_insertion.bam"
samtools index "$tmp/leftshifted_insertion.bam"
"$bin" \
  -q \
  -r "$tmp/poly.fa" \
  -s S1 \
  --phase-from-bam "$tmp/leftshifted_insertion.bam" \
  --phase-realign-overhang 1 \
  --emit all-sites \
  "$tmp/leftshifted_insertion.vcf" > "$tmp/leftshifted_insertion.out.vcf"
grep -v '^#' "$tmp/leftshifted_insertion.out.vcf" > "$tmp/leftshifted_insertion.body.vcf"
cat > "$tmp/leftshifted_insertion.expected.vcf" <<'VCF'
chrR	2	.	A	C	.	PASS	.	GT:PS	0|1:2
chrR	20	.	A	AA	.	PASS	.	GT:PS	0|1:2
VCF
diff -u "$tmp/leftshifted_insertion.expected.vcf" "$tmp/leftshifted_insertion.body.vcf"

else
  echo "samtools not found; skipping generated BAM phasing regression tests"
fi

"$bin" \
  -r "$ref" \
  -s S1 \
  --phase-from-bam "$fixtures/read_phase.bam" \
  --no-header \
  "$fixtures/read_phase.vcf" > "$tmp/out.with-summary.vcf" 2> "$tmp/summary.err"

grep -q 'phase_mnv: bam_phase input=' "$tmp/summary.err"
grep -q 'algorithm=mec' "$tmp/summary.err"
grep -q 'selected_reads=4' "$tmp/summary.err"
grep -q 'candidates=4' "$tmp/summary.err"
grep -q 'phased_variants=4' "$tmp/summary.err"
grep -q 'emitted_calls=2' "$tmp/summary.err"

echo "phase_mnv BAM phasing tests passed"
