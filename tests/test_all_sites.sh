#!/usr/bin/env bash
set -euo pipefail

bin=${1:-./target/release/phase_mnv_rs}
root=$(cd "$(dirname "$0")/.." && pwd)
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

"$bin" \
  --emit all-sites \
  --threads 2 \
  --reference "$root/tests/fixtures/ref.fa" \
  --phase-from-bam "$root/tests/fixtures/read_phase.bam" \
  "$root/tests/fixtures/read_phase.vcf" \
  > "$tmp/all_sites.vcf" \
  2> "$tmp/all_sites.stderr"

grep -q '^##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">$' "$tmp/all_sites.vcf"
grep -q '^##phase_mnv_header_policy=preserve_input_header_and_append_phase_mnv_records$' "$tmp/all_sites.vcf"
grep -q '^##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set assigned by phase_mnv_rs">$' "$tmp/all_sites.vcf"
grep -q 'output_format=vcf threads=2 emit=all-sites' "$tmp/all_sites.stderr"

cat > "$tmp/expected.body" <<'EOF'
chr1	1	.	A	G	.	PASS	.	GT:PS	0|1:1
chr1	2	.	C	T	.	PASS	.	GT:PS	0|1:1
chr1	4	.	T	C	.	PASS	.	GT:PS	1|0:1
chr1	5	.	A	G	.	PASS	.	GT:PS	1|0:1
EOF

grep -v '^#' "$tmp/all_sites.vcf" > "$tmp/all_sites.body"
diff -u "$tmp/expected.body" "$tmp/all_sites.body"

if command -v bcftools >/dev/null 2>&1; then
  "$bin" \
    --emit all-sites \
    --quiet \
    --threads 2 \
    --reference "$root/tests/fixtures/ref.fa" \
    --phase-from-bam "$root/tests/fixtures/read_phase.bam" \
    --output "$tmp/all_sites.bcf" \
    "$root/tests/fixtures/read_phase.vcf"
  bcftools view -h "$tmp/all_sites.bcf" > "$tmp/all_sites.bcf.header"
  grep -q '^##phase_mnv_header_policy=preserve_input_header_and_append_phase_mnv_records$' "$tmp/all_sites.bcf.header"
  bcftools view -H "$tmp/all_sites.bcf" > "$tmp/all_sites.bcf.body"
  diff -u "$tmp/expected.body" "$tmp/all_sites.bcf.body"
fi
