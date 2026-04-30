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
  --region chr1:1-6 \
  --position-tsv "$tmp/positions.tsv" \
  "$fixtures/read_phase.bam" > "$tmp/model.with-positions.tsv"
grep -qx $'read_pos\tbaseq_group\tobservations\tmatches\tmismatches\tinsertions\tdeletions\terror_rate\tempirical_q' "$tmp/positions.tsv"
grep -qx $'1\tall\t4\t2\t2\t0\t0\t0.500000\t3.010' "$tmp/positions.tsv"
grep -qx $'1\thigh\t4\t2\t2\t0\t0\t0.500000\t3.010' "$tmp/positions.tsv"
grep -qx $'3\tall\t4\t4\t0\t0\t0\t0.000000\tinf' "$tmp/positions.tsv"

grep -qx $'overall\tall\t24\t16\t8\t0\t0\t0.333333\t4.771' "$tmp/model.with-positions.tsv"

cat > "$tmp/rev.fa" <<'EOF'
>chr1
ACGT
EOF
printf 'chr1\t4\t6\t4\t5\n' > "$tmp/rev.fa.fai"
cat > "$tmp/rev.sam" <<'EOF'
@HD	VN:1.6	SO:unknown
@SQ	SN:chr1	LN:4
rev	16	chr1	1	60	4M	*	0	0	TCGT	IIII
EOF
"$bin" --reference "$tmp/rev.fa" --position-tsv "$tmp/rev.positions.tsv" "$tmp/rev.sam" > "$tmp/rev.tsv"
grep -qx $'1\tall\t1\t1\t0\t0\t0\t0.000000\tinf' "$tmp/rev.positions.tsv"
grep -qx $'4\tall\t1\t0\t1\t0\t0\t1.000000\t-0.000' "$tmp/rev.positions.tsv"

cat > "$tmp/indel.sam" <<'EOF'
@HD	VN:1.6	SO:unknown
@SQ	SN:chr1	LN:4
del	0	chr1	1	60	1M1D2M	*	0	0	AGT	III
EOF
"$bin" --reference "$tmp/rev.fa" --position-tsv "$tmp/indel.positions.tsv" "$tmp/indel.sam" > "$tmp/indel.tsv"
grep -qx $'overall\tall\t4\t3\t0\t0\t1\t0.250000\t6.021' "$tmp/indel.tsv"
grep -qx $'baseq\t40-49\t3\t3\t0\t0\t0\t0.000000\tinf' "$tmp/indel.tsv"
grep -qx $'1\tall\t2\t1\t0\t0\t1\t0.500000\t3.010' "$tmp/indel.positions.tsv"

if [[ -e /dev/full ]]; then
  if "$bin" --reference "$fixtures/ref.fa" --position-tsv /dev/full "$fixtures/read_phase.bam" > "$tmp/full.out" 2> "$tmp/full.err"; then
    echo "bam_error_model unexpectedly succeeded when --position-tsv could not be written" >&2
    exit 1
  fi
  grep -q 'failed to .*position TSV' "$tmp/full.err"
fi

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
