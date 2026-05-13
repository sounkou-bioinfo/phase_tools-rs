#!/usr/bin/env bash
set -euo pipefail

bin=${1:?usage: tests/test_negative.sh /path/to/phase_mnv_binary}
repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
fixtures="$repo_root/tests/fixtures"
ref="$fixtures/ref.fa"
valid_vcf="$fixtures/phased_mnv.vcf"
truncated_vcf="$fixtures/truncated.vcf.gz"

if [[ ! -x "$bin" ]]; then
  echo "error: binary is not executable: $bin" >&2
  exit 2
fi

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

expect_fail() {
  local name=$1
  local pattern=$2
  shift 2
  local out="$tmp/${name//[^A-Za-z0-9_.-]/_}.out"
  local err="$tmp/${name//[^A-Za-z0-9_.-]/_}.err"
  set +e
  "$@" > "$out" 2> "$err"
  local status=$?
  set -e
  if [[ "$status" -eq 0 ]]; then
    echo "negative test failed: $name unexpectedly exited 0" >&2
    echo "command: $*" >&2
    echo "stdout:" >&2
    sed -n '1,80p' "$out" >&2
    echo "stderr:" >&2
    sed -n '1,80p' "$err" >&2
    exit 1
  fi
  if ! grep -Eiq "$pattern" "$err"; then
    echo "negative test failed: $name stderr did not match /$pattern/" >&2
    echo "command: $*" >&2
    echo "exit status: $status" >&2
    echo "stdout:" >&2
    sed -n '1,80p' "$out" >&2
    echo "stderr:" >&2
    sed -n '1,120p' "$err" >&2
    exit 1
  fi
}

missing_input="$tmp/does-not-exist.vcf"
missing_ref="$tmp/does-not-exist.fa"
for required_fixture in "$ref" "$valid_vcf" "$fixtures/ref_mismatch.vcf" "$truncated_vcf"; do
  if [[ ! -f "$required_fixture" ]]; then
    echo "error: missing negative-test fixture: $required_fixture" >&2
    exit 2
  fi
done

expect_fail "missing required reference" \
  "reference is required" \
  "$bin" "$valid_vcf"

expect_fail "missing input VCF" \
  "cannot open input|failed to read VCF/BCF header|No such file|Failed to open" \
  "$bin" -r "$ref" -s S1 "$missing_input"

expect_fail "missing FASTA reference" \
  "FASTA|reference|fai|index|No such file|cannot load" \
  "$bin" -r "$missing_ref" -s S1 "$valid_vcf"

expect_fail "missing sample" \
  "sample .*not found|Available samples" \
  "$bin" -r "$ref" -s DOES_NOT_EXIST "$valid_vcf"

expect_fail "negative max gap" \
  "max-gap must be >= 0" \
  "$bin" -r "$ref" --max-gap -1 -s S1 "$valid_vcf"

expect_fail "write-index plain VCF" \
  "write-index cannot index plain VCF" \
  "$bin" -r "$ref" -s S1 --write-index -o "$tmp/out.vcf" "$valid_vcf"

expect_fail "write-index tbi bcf" \
  "BCF output requires CSI" \
  "$bin" -r "$ref" -s S1 --write-index=tbi -o "$tmp/out.bcf" "$valid_vcf"

expect_fail "combined no-header" \
  "combined preserves the original VCF/BCF header" \
  "$bin" -r "$ref" -s S1 --emit combined --no-header "$valid_vcf"

expect_fail "combined phase-from-bam" \
  "combined currently requires input-phased VCF/BCF" \
  "$bin" -r "$ref" -s S1 --emit combined --phase-from-bam "$fixtures/read_phase.bam" "$valid_vcf"

expect_fail "REF FASTA mismatch" \
  "REF/FASTA mismatch|Use --no-ref-check" \
  "$bin" -r "$ref" -s S1 "$fixtures/ref_mismatch.vcf"

expect_fail "unsupported selected allele fail policy" \
  "unsupported selected ALT allele|kind=symbolic_or_breakend" \
  "$bin" -r "$ref" -s S1 --max-gap 1 --unsupported-alleles fail "$fixtures/symbolic.vcf"

expect_fail "truncated gzip VCF" \
  "truncated|failed to read VCF/BCF header|failed to read VCF/BCF record|Reading GZIP stream failed" \
  "$bin" -r "$ref" -s S1 "$truncated_vcf"

echo "phase_mnv negative tests passed for $bin"
