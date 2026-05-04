#!/usr/bin/env bash
set -euo pipefail

bin=${1:?usage: $0 <bam_ancestry_binary>}
repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
fixtures="$repo_root/tests/fixtures"
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

cat > "$tmp/anchors.tsv" <<'EOF'
chrom	pos	ref	alt	POP1	POP2
chr1	1	A	G	0.0	1.0
chr1	2	C	T	0.0	1.0
chr1	4	T	C	0.0	1.0
EOF

"$bin" \
  --reference "$fixtures/ref.fa" \
  --bam "$fixtures/read_phase.bam" \
  --anchors "$tmp/anchors.tsv" > "$tmp/out.tsv"

grep -qx $'#anchors\t3' "$tmp/out.tsv"
grep -qx $'#used_anchors\t3' "$tmp/out.tsv"
grep -qx $'#weighted_sse\t0.000000' "$tmp/out.tsv"
grep -qx $'population\tproportion' "$tmp/out.tsv"
grep -qx $'POP1\t0.500000' "$tmp/out.tsv"
grep -qx $'POP2\t0.500000' "$tmp/out.tsv"
grep -qx $'chr1\t1\tA\tG\t1\t4\t2\t2\t0\t0\t0.500000\t0.500000\t0.000000\t0.000000\t1.000000' "$tmp/out.tsv"

"$bin" \
  --reference "$fixtures/ref.fa" \
  --bam "$fixtures/read_phase.bam" \
  --anchors "$tmp/anchors.tsv" \
  --populations POP2,POP1 > "$tmp/reordered.tsv"
grep -qx $'POP2\t0.500000' "$tmp/reordered.tsv"
grep -qx $'POP1\t0.500000' "$tmp/reordered.tsv"

echo -e 'chrom\tpos\tref\talt\tPOP1\tPOP2\nchr1\t1\tC\tG\t0.0\t1.0' > "$tmp/ref_mismatch.tsv"
if "$bin" --reference "$fixtures/ref.fa" --bam "$fixtures/read_phase.bam" --anchors "$tmp/ref_mismatch.tsv" > "$tmp/ref_mismatch.out" 2> "$tmp/ref_mismatch.err"; then
  echo "bam_ancestry unexpectedly accepted anchor/FASTA REF mismatch" >&2
  exit 1
fi
grep -q 'anchor REF mismatch' "$tmp/ref_mismatch.err"

"$bin" \
  --reference "$fixtures/ref.fa" \
  --bam "$fixtures/read_phase.bam" \
  --anchors "$tmp/anchors.tsv" \
  --min-baseq 41 > "$tmp/baseq.tsv" 2> "$tmp/baseq.err" && {
    echo "bam_ancestry unexpectedly fit anchors with no observations" >&2
    exit 1
  }
grep -q 'no anchors had enough REF+ALT observations' "$tmp/baseq.err"

echo "bam_ancestry tests passed"
