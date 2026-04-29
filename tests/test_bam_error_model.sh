#!/usr/bin/env bash
set -euo pipefail

bin=${1:?usage: $0 <bam_error_model_binary>}
repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
fixtures="$repo_root/tests/fixtures"
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

"$bin" \
  --reference "$fixtures/ref.fa" \
  --region chr1:1-6 \
  --region chr1:1-6 \
  --max-reads 10 \
  "$fixtures/read_phase.bam" > "$tmp/model.tsv"

grep -qx $'#reads\t4' "$tmp/model.tsv"
grep -qx $'overall\tall\t24\t16\t8\t0\t0\t0.333333\t4.771' "$tmp/model.tsv"
grep -qx $'baseq\t40-49\t24\t16\t8\t0\t0\t0.333333\t4.771' "$tmp/model.tsv"
grep -qx $'mapq\t60+\t24\t16\t8\t0\t0\t0.333333\t4.771' "$tmp/model.tsv"

"$bin" \
  --reference "$fixtures/ref.fa" \
  --region chr1:1-1 \
  "$fixtures/read_phase.bam" > "$tmp/model.clipped.tsv"
grep -qx $'#reads\t4' "$tmp/model.clipped.tsv"
grep -qx $'overall\tall\t4\t2\t2\t0\t0\t0.500000\t3.010' "$tmp/model.clipped.tsv"

"$bin" \
  --reference "$fixtures/ref.fa" \
  --region chr1:1-1 \
  --region chr1:6-6 \
  "$fixtures/read_phase.bam" > "$tmp/model.disjoint.tsv"
grep -qx $'#reads\t4' "$tmp/model.disjoint.tsv"
grep -qx $'overall\tall\t8\t6\t2\t0\t0\t0.250000\t6.021' "$tmp/model.disjoint.tsv"

"$bin" \
  --reference "$fixtures/ref.fa" \
  --region chr1:1-1 \
  --region chr1:6-6 \
  --max-reads 2 \
  "$fixtures/read_phase.bam" > "$tmp/model.disjoint.max.tsv"
grep -qx $'#reads\t2' "$tmp/model.disjoint.max.tsv"
grep -qx $'overall\tall\t4\t2\t2\t0\t0\t0.500000\t3.010' "$tmp/model.disjoint.max.tsv"

"$bin" \
  --reference "$fixtures/ref.fa" \
  --region chr1:1-6 \
  --min-mapq 61 \
  "$fixtures/read_phase.bam" > "$tmp/model.mapq.tsv"
grep -qx $'#reads\t0' "$tmp/model.mapq.tsv"
grep -qx $'overall\tall\t0\t0\t0\t0\t0\tNA\tNA' "$tmp/model.mapq.tsv"

cat > "$tmp/mapq255.sam" <<'EOF'
@HD	VN:1.6	SO:unknown
@SQ	SN:chr1	LN:12
mapq_unknown	0	chr1	1	255	4M	*	0	0	ACGT	IIII
EOF
"$bin" --reference "$fixtures/ref.fa" --max-reads 1 "$tmp/mapq255.sam" > "$tmp/mapq255.tsv"
grep -qx $'#reads\t1' "$tmp/mapq255.tsv"
grep -qx $'mapq\tunknown\t4\t4\t0\t0\t0\t0.000000\tinf' "$tmp/mapq255.tsv"
"$bin" --reference "$fixtures/ref.fa" --min-mapq 1 "$tmp/mapq255.sam" > "$tmp/mapq255.filtered.tsv"
grep -qx $'#reads\t0' "$tmp/mapq255.filtered.tsv"

echo "bam_error_model tests passed"
