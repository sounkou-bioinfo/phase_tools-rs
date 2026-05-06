#!/usr/bin/env bash
set -euo pipefail

bin=${1:?usage: $0 <phase_adjudicate_binary>}
repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
fixtures="$repo_root/tests/fixtures"
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

cat > "$tmp/pairs.tsv" <<'EOF'
chrom	prev_pos	pos	truth_ps	query_ps	prev_orientation	orientation	status	prev_gt_truth	gt_truth	prev_gt_query	gt_query
chr1	1	2	1	1	0	1	switch	0|1	0|1	0|1	1|0
chr1	2	4	1	1	0	1	switch	0|1	1|0	0|1	0|1
chr1	4	5	1	1	0	0	match	1|0	1|0	0|1	0|1
EOF

"$bin" \
  --reference "$fixtures/ref.fa" \
  --bam "$fixtures/read_phase.bam" \
  --variants "$fixtures/read_phase.vcf" \
  --pair-tsv "$tmp/pairs.tsv" > "$tmp/out.tsv"

grep -qx $'chrom\tprev_pos\tpos\ttruth_parity\tquery_parity\tusable_reads\tspanning_reads\tinformative_reads\ttruth_support\tquery_support\tother_support\tmean_min_baseq\tmean_mapq\tmapq_unknown_reads\tforward_reads\treverse_reads\twinner\tambiguous\treason' "$tmp/out.tsv"
grep -qx $'chr1\t1\t2\t0\t1\t4\t4\t4\t4\t0\t0\t40.000\t60.000\t0\t4\t0\ttruth\tfalse\tevidence' "$tmp/out.tsv"
grep -qx $'chr1\t2\t4\t1\t0\t4\t4\t4\t4\t0\t0\t40.000\t60.000\t0\t4\t0\ttruth\tfalse\tevidence' "$tmp/out.tsv"
grep -qx $'chr1\t4\t5\t0\t0\t4\t4\t4\t4\t4\t0\t40.000\t60.000\t0\t4\t0\tboth\tfalse\tsame_phase' "$tmp/out.tsv"

case "$(uname -m)" in
  x86_64|amd64)
    if ! command -v samtools >/dev/null 2>&1; then
      echo "samtools not found; skipping phase_adjudicate assembly sidecar test"
    else
    python3 - <<'PY' > "$tmp/asm_ref.fa"
ref = "ACGTGCTAGCTAGGATCCGATCGATGCTAGCTAGCATCGATCGTTAGCTAGGCTAACCGGTTAACCGGTTAGCATCGATCG"
print(">chrA")
print(ref)
PY
    samtools faidx "$tmp/asm_ref.fa"
    python3 - <<'PY' > "$tmp/asm.sam"
ref = "ACGTGCTAGCTAGGATCCGATCGATGCTAGCTAGCATCGATCGTTAGCTAGGCTAACCGGTTAACCGGTTAGCATCGATCG"
hap = list(ref)
hap[9] = "T"
hap[24] = "A"
hap = "".join(hap)
print("@HD\tVN:1.6\tSO:coordinate")
print(f"@SQ\tSN:chrA\tLN:{len(ref)}")
for idx, start in enumerate([0, 8, 16, 24, 32, 40], 1):
    seq = hap[start:start + 40]
    print(f"asm{idx}\t0\tchrA\t{start + 1}\t60\t{len(seq)}M\t*\t0\t0\t{seq}\t" + "I" * len(seq))
PY
    samtools view -b "$tmp/asm.sam" | samtools sort -o "$tmp/asm.bam"
    samtools index "$tmp/asm.bam"
    cat > "$tmp/asm.vcf" <<'EOF'
##fileformat=VCFv4.3
##contig=<ID=chrA>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
chrA	10	.	C	T	.	PASS	.	GT	0/1
chrA	25	.	T	A	.	PASS	.	GT	0/1
EOF
    cat > "$tmp/asm_pairs.tsv" <<'EOF'
chrom	prev_pos	pos	prev_gt_truth	gt_truth	prev_gt_query	gt_query
chrA	10	25	0|1	0|1	0|1	1|0
EOF
    "$bin" \
      --reference "$tmp/asm_ref.fa" \
      --bam "$tmp/asm.bam" \
      --variants "$tmp/asm.vcf" \
      --pair-tsv "$tmp/asm_pairs.tsv" \
      --assembly-fasta "$tmp/assembly.fa" \
      --assembly-tsv "$tmp/assembly.tsv" \
      --assembly-window 30 \
      --assembly-context 25 \
      --assembly-min-asm-ovlp 12 > "$tmp/out.assembly.tsv"
    grep -q '^>chrA:10-25|unitig=1' "$tmp/assembly.fa"
    grep -qx $'chrom\tprev_pos\tpos\tunitig\tinput_reads\tunitig_len\tsupporting_reads\tassembly_start\tassembly_end\tbest_prev_allele\tbest_allele\tbest_parity\tbest_distance\tsecond_distance\tstatus\tsupports_truth\tsupports_query' "$tmp/assembly.tsv"
    grep -Eq $'^chrA\t10\t25\t1\t6\t[0-9]+\t[0-9]+\t1\t50\t1\t1\t0\t0\t[0-9]+\tinformative\ttrue\tfalse$' "$tmp/assembly.tsv"
    "$bin" \
      --reference "$tmp/asm_ref.fa" \
      --bam "$tmp/asm.bam" \
      --variants "$tmp/asm.vcf" \
      --pair-tsv "$tmp/asm_pairs.tsv" \
      --min-baseq 41 \
      --assembly-tsv "$tmp/assembly.decision.tsv" \
      --use-assembly-decision \
      --assembly-window 30 \
      --assembly-context 25 \
      --assembly-min-asm-ovlp 12 > "$tmp/out.assembly_decision.tsv"
    grep -qx $'chrA\t10\t25\t0\t1\t4\t0\t0\t0\t0\t0\tNA\tNA\t0\t0\t0\ttruth\tfalse\tassembly_evidence' "$tmp/out.assembly_decision.tsv"
    fi
    ;;
  *) echo "phase_adjudicate assembly sidecar test skipped on non-x86_64 host" ;;
esac

"$bin" \
  --reference "$fixtures/ref.fa" \
  --bam "$fixtures/read_phase.bam" \
  --variants "$fixtures/read_phase.vcf" \
  --pair-tsv "$tmp/pairs.tsv" \
  --min-baseq 41 > "$tmp/baseq.tsv"
grep -qx $'chr1\t1\t2\t0\t1\t4\t0\t0\t0\t0\t0\tNA\tNA\t0\t0\t0\tnone\ttrue\tno_informative_reads' "$tmp/baseq.tsv"
grep -qx $'chr1\t4\t5\t0\t0\t4\t0\t0\t0\t0\t0\tNA\tNA\t0\t0\t0\tboth\tfalse\tsame_phase' "$tmp/baseq.tsv"

cat > "$tmp/unsupported_pairs.tsv" <<'EOF'
chrom	prev_pos	pos	prev_gt_truth	gt_truth	prev_gt_query	gt_query
chr1	1	2	0/1	0|1	0|1	0|1
chr1	1	2	0|1	0|1	0|0	0|1
EOF
"$bin" \
  --reference "$fixtures/ref.fa" \
  --bam "$fixtures/read_phase.bam" \
  --variants "$fixtures/read_phase.vcf" \
  --pair-tsv "$tmp/unsupported_pairs.tsv" > "$tmp/unsupported.tsv"
grep -qx $'chr1\t1\t2\tNA\tNA\t0\t0\t0\t0\t0\t0\tNA\tNA\t0\t0\t0\tnone\ttrue\tunsupported_truth_gt' "$tmp/unsupported.tsv"
grep -qx $'chr1\t1\t2\t0\tNA\t0\t0\t0\t0\t0\t0\tNA\tNA\t0\t0\t0\tnone\ttrue\tunsupported_query_gt' "$tmp/unsupported.tsv"

cat > "$tmp/bad_pairs.tsv" <<'EOF'
chrom	prev_pos	pos	prev_gt_truth	gt_truth	prev_gt_query	gt_query
chr1	1	3	0|1	0|1	0|1	0|1
EOF
"$bin" \
  --reference "$fixtures/ref.fa" \
  --bam "$fixtures/read_phase.bam" \
  --variants "$fixtures/read_phase.vcf" \
  --pair-tsv "$tmp/bad_pairs.tsv" > "$tmp/missing.tsv"
grep -qx $'chr1\t1\t3\t0\t0\t0\t0\t0\t0\t0\t0\tNA\tNA\t0\t0\t0\tnone\ttrue\tmissing_variant' "$tmp/missing.tsv"

cat > "$tmp/duplicate_pos.vcf" <<'EOF'
##fileformat=VCFv4.3
##contig=<ID=chr1>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	1	.	A	G	.	PASS	.
chr1	1	.	A	C	.	PASS	.
chr1	2	.	C	T	.	PASS	.
EOF
if "$bin" \
  --reference "$fixtures/ref.fa" \
  --bam "$fixtures/read_phase.bam" \
  --variants "$tmp/duplicate_pos.vcf" \
  --pair-tsv "$tmp/pairs.tsv" > "$tmp/duplicate.out" 2> "$tmp/duplicate.err"; then
  echo "phase_adjudicate unexpectedly accepted duplicate-position SNVs" >&2
  exit 1
fi
grep -q 'duplicate biallelic SNV records at chr1:1' "$tmp/duplicate.err"

echo "phase_adjudicate tests passed"
