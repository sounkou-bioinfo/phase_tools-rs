#!/usr/bin/env bash
set -euo pipefail

# Fast WhatsHap-vs-Rust phasing comparison without hap.py.
#
# truth path: input VCF -> unphase -> external `whatshap phase`
# query path: input VCF -> `phase_mnv_rs --emit all-sites --phase-from-bam`
# comparison: `phase_compare` phase-block/pair/switch statistics

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)

phase_mnv_bin=${PHASE_MNV_BIN:-"$repo_root/target/release/phase_mnv_rs"}
phase_compare_bin=${PHASE_COMPARE_BIN:-"$repo_root/target/release/phase_compare"}
ref=${REF:-"$repo_root/tests/fixtures/ref.fa"}
vcf=${VCF:-"$repo_root/tests/fixtures/read_phase.vcf"}
bam=${BAM:-"$repo_root/tests/fixtures/read_phase.bam"}
sample=${SAMPLE:-S1}
phase_algorithm=${PHASE_ALGORITHM:-mec}
phase_max_coverage=${PHASE_MAX_COVERAGE:-15}
whatshap_env=${WHATSHAP_ENV:-phase-mnv-whatshap}
threads=${THREADS:-1}
allow_nonperfect=${ALLOW_NONPERFECT:-0}
max_switch_errors=${MAX_SWITCH_ERRORS:-0}
max_switch_rate=${MAX_SWITCH_RATE:-0}
keep_tmp=${KEEP_TMP:-0}

usage() {
  cat <<'USAGE'
usage: scripts/compare_whatshap_phase.sh

Compares external WhatsHap phasing against Rust phase_mnv_rs BAM phasing using
the fast native phase_compare binary. No hap.py is used.

Environment overrides:
  PHASE_MNV_BIN       phase_mnv_rs binary (default: target/release/phase_mnv_rs)
  PHASE_COMPARE_BIN   phase_compare binary (default: target/release/phase_compare)
  REF                 FASTA reference (default: tests/fixtures/ref.fa)
  VCF                 input VCF/VCF.GZ/BCF (default: tests/fixtures/read_phase.vcf)
  BAM                 indexed BAM/CRAM (default: tests/fixtures/read_phase.bam)
  SAMPLE              sample name (default: S1)
  PHASE_ALGORITHM     Rust --phase-algorithm for query path (default: mec)
  PHASE_MAX_COVERAGE  Rust --phase-max-coverage for query path (default: 15)
  WHATSHAP_BIN        explicit whatshap command/path
  WHATSHAP_ENV        micromamba env fallback for whatshap (default: phase-mnv-whatshap)
  THREADS             htslib reader threads for phase_compare/Rust IO (default: 1)
  MAX_SWITCH_ERRORS   accepted switch-error count (default: 0)
  MAX_SWITCH_RATE     accepted switch-error rate (default: 0)
  ALLOW_NONPERFECT=1  report summary but do not fail thresholds
  KEEP_TMP=1          keep temporary comparison directory

Generated comparison VCFs are sanitised with scripts/sanitize_vcf_headers.py so
bcftools/producer command headers and local path-bearing phase_mnv/reference
records are removed before phase_compare.
USAGE
}

if [[ ${1:-} == "-h" || ${1:-} == "--help" ]]; then
  usage
  exit 0
fi

die() {
  echo "error: $*" >&2
  exit 1
}

require_file() {
  [[ -e "$1" ]] || die "missing required file: $1"
}

make_cmd() {
  local explicit=$1
  local binary=$2
  local env_name=$3
  local -n out=$4
  local env_var=$5
  if [[ -n "$explicit" ]]; then
    # Intentional shell splitting permits wrappers such as "micromamba run -n env whatshap".
    # shellcheck disable=SC2206
    out=($explicit)
    return
  fi
  if command -v "$binary" >/dev/null 2>&1; then
    out=("$binary")
    return
  fi
  if command -v micromamba >/dev/null 2>&1 && micromamba run -n "$env_name" "$binary" --help >/dev/null 2>&1; then
    out=(micromamba run -n "$env_name" "$binary")
    return
  fi
  die "$binary not found; set $env_var or create micromamba env '$env_name'"
}

require_file "$phase_mnv_bin"
require_file "$phase_compare_bin"
require_file "$ref"
require_file "$vcf"
require_file "$bam"

whatshap_cmd=()
make_cmd "${WHATSHAP_BIN:-}" whatshap "$whatshap_env" whatshap_cmd WHATSHAP_BIN

tmp=${TMPDIR:-/tmp}/phase-mnv-whatshap-phase.$$
mkdir -p "$tmp"
if [[ "$keep_tmp" != "1" ]]; then
  trap 'rm -rf "$tmp"' EXIT
else
  echo "compare_whatshap_phase: keeping tmp=$tmp" >&2
fi

unphased_vcf="$tmp/input.unphased.vcf"
whatshap_raw="$tmp/whatshap.raw.vcf"
whatshap_vcf="$tmp/whatshap.sanitized.vcf"
rust_raw="$tmp/rust.all_sites.raw.vcf"
rust_vcf="$tmp/rust.all_sites.sanitized.vcf"
summary="$tmp/phase_compare.summary.tsv"
switch_bed="$tmp/phase_compare.switches.bed"
pair_tsv="$tmp/phase_compare.pairs.tsv"

python3 "$repo_root/scripts/unphase_vcf.py" "$vcf" > "$unphased_vcf"

"${whatshap_cmd[@]}" phase \
  --reference "$ref" \
  --sample "$sample" \
  -o "$whatshap_raw" \
  "$unphased_vcf" \
  "$bam" > "$tmp/whatshap.log" 2>&1
python3 "$repo_root/scripts/sanitize_vcf_headers.py" "$whatshap_raw" "$whatshap_vcf"

"$phase_mnv_bin" \
  -q \
  --emit all-sites \
  --reference "$ref" \
  --sample "$sample" \
  --phase-from-bam "$bam" \
  --phase-algorithm "$phase_algorithm" \
  --phase-max-coverage "$phase_max_coverage" \
  --threads "$threads" \
  --output "$rust_raw" \
  "$vcf" > "$tmp/rust.all_sites.log" 2>&1
python3 "$repo_root/scripts/sanitize_vcf_headers.py" "$rust_raw" "$rust_vcf"

"$phase_compare_bin" \
  --sample "$sample" \
  --threads "$threads" \
  --switch-bed "$switch_bed" \
  --pair-tsv "$pair_tsv" \
  "$whatshap_vcf" \
  "$rust_vcf" > "$summary"

cat "$summary"

total=$(awk -F'\t' '$1=="TOTAL" {print}' "$summary")
[[ -n "$total" ]] || die "phase_compare TOTAL row missing"
switch_errors=$(awk -F'\t' '$1=="TOTAL" {print $16}' "$summary")
switch_rate=$(awk -F'\t' '$1=="TOTAL" {print $17}' "$summary")

if awk -v e="$switch_errors" -v maxe="$max_switch_errors" -v r="$switch_rate" -v maxr="$max_switch_rate" 'BEGIN { exit !((e+0) <= (maxe+0) && (r == "NA" || (r+0) <= (maxr+0))) }'; then
  echo "compare_whatshap_phase: phase comparison passed"
else
  echo "compare_whatshap_phase: switch thresholds exceeded: switch_errors=$switch_errors switch_rate=$switch_rate" >&2
  echo "temporary files: $tmp" >&2
  if [[ "$allow_nonperfect" == "1" ]]; then
    exit 0
  fi
  exit 1
fi
