#!/usr/bin/env bash
set -euo pipefail

# Compare phase_mnv_rs against vcflib's vcfgeno2haplo on a narrow, fair slice:
# adjacent phased SNV blocks in one sample.  This is a semantic/projection
# comparison, not a byte-identical comparison, because vcfgeno2haplo emits
# haplotype-allele VCF records and may pass through non-cluster input records,
# while phase_mnv_rs emits normalized TYPE=MNV/TYPE=COMPLEX records.

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
rs_bin=${RS_BIN:-"$repo_root/target/release/phase_mnv_rs"}
ref=${REF:-"$repo_root/tests/fixtures/ref.fa"}
vcf=${VCF:-"$repo_root/tests/fixtures/vcfgeno2haplo_compare.vcf"}
sample=${SAMPLE:-S1}
window=${WINDOW:-1}
env_name=${VCFLIB_ENV:-phase-mnv-vcflib}
auto_create=${VCFLIB_AUTO_CREATE:-1}
vcflib_spec=${VCFLIB_SPEC:-vcflib=1.0.15}
export MAMBA_ROOT_PREFIX=${MAMBA_ROOT_PREFIX:-$HOME/micromamba}

usage() {
  cat <<'USAGE'
usage: scripts/compare_vcfgeno2haplo.sh

Environment variables:
  RS_BIN                 phase_mnv_rs binary (default: target/release/phase_mnv_rs)
  VCF                    input VCF fixture (default: tests/fixtures/vcfgeno2haplo_compare.vcf)
  REF                    FASTA reference (default: tests/fixtures/ref.fa)
  SAMPLE                 sample name (default: S1)
  WINDOW                 vcfgeno2haplo window size (default: 1)
  VCFGENO2HAPLO_BIN      explicit vcfgeno2haplo executable
  VCFLIB_ENV             micromamba env name (default: phase-mnv-vcflib)
  VCFLIB_SPEC            micromamba package spec (default: vcflib=1.0.15)
  VCFLIB_AUTO_CREATE     create micromamba env if needed (default: 1)
USAGE
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

if [[ ! -x "$rs_bin" ]]; then
  echo "error: phase_mnv_rs binary is not executable: $rs_bin" >&2
  echo "hint: run 'make release' first" >&2
  exit 2
fi
if [[ ! -f "$vcf" ]]; then
  echo "error: input VCF does not exist: $vcf" >&2
  exit 2
fi
if [[ ! -f "$ref" ]]; then
  echo "error: reference FASTA does not exist: $ref" >&2
  exit 2
fi

vcfgeno2haplo_cmd=()
ensure_vcfgeno2haplo() {
  if [[ -n "${VCFGENO2HAPLO_BIN:-}" ]]; then
    if [[ ! -x "$VCFGENO2HAPLO_BIN" ]]; then
      echo "error: VCFGENO2HAPLO_BIN is not executable: $VCFGENO2HAPLO_BIN" >&2
      exit 2
    fi
    vcfgeno2haplo_cmd=("$VCFGENO2HAPLO_BIN")
    return
  fi

  if command -v vcfgeno2haplo >/dev/null 2>&1; then
    vcfgeno2haplo_cmd=(vcfgeno2haplo)
    return
  fi

  if ! command -v micromamba >/dev/null 2>&1; then
    echo "error: vcfgeno2haplo not found and micromamba is not installed" >&2
    exit 2
  fi

  if ! micromamba run -n "$env_name" vcfgeno2haplo --help >/dev/null 2>&1; then
    if [[ "$auto_create" != "1" ]]; then
      echo "error: micromamba env '$env_name' does not provide vcfgeno2haplo" >&2
      echo "hint: set VCFLIB_AUTO_CREATE=1 or install vcflib into that env" >&2
      exit 2
    fi
    echo "Creating micromamba env '$env_name' with $vcflib_spec..." >&2
    micromamba create -y -n "$env_name" -c conda-forge -c bioconda "$vcflib_spec" >&2
  fi
  vcfgeno2haplo_cmd=(micromamba run -n "$env_name" vcfgeno2haplo)
}

project_phase_mnv() {
  awk -F'\t' 'BEGIN { OFS="\t" }
    /^#/ { next }
    {
      split($10, sample_fields, ":")
      print $1, $2, $4, $5, sample_fields[1]
    }' "$1"
}

project_vcfgeno2haplo() {
  awk -F'\t' 'BEGIN { OFS="\t" }
    /^#/ { next }
    length($4) > 1 || length($5) > 1 {
      split($10, sample_fields, ":")
      print $1, $2, $4, $5, sample_fields[1]
    }' "$1"
}

ensure_vcfgeno2haplo

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

phase_out="$tmp/phase_mnv.vcf"
vcflib_out="$tmp/vcfgeno2haplo.vcf"
phase_projection="$tmp/phase_mnv.projected.tsv"
vcflib_projection="$tmp/vcfgeno2haplo.projected.tsv"

"$rs_bin" -q -r "$ref" -s "$sample" --no-header -o "$phase_out" "$vcf"
"${vcfgeno2haplo_cmd[@]}" -r "$ref" -w "$window" "$vcf" > "$vcflib_out"

project_phase_mnv "$phase_out" > "$phase_projection"
project_vcfgeno2haplo "$vcflib_out" > "$vcflib_projection"

if ! diff -u "$vcflib_projection" "$phase_projection"; then
  cat >&2 <<EOF
vcfgeno2haplo semantic projection differs from phase_mnv_rs.
Compared columns: CHROM POS REF ALT GT.
vcfgeno2haplo output: $vcflib_out
phase_mnv_rs output: $phase_out
EOF
  exit 1
fi

vcflib_version=$("${vcfgeno2haplo_cmd[@]}" --version 2>&1 | head -1 || true)
if [[ -z "$vcflib_version" ]]; then
  vcflib_version="vcfgeno2haplo version unavailable"
fi

cat <<EOF
vcfgeno2haplo comparison passed
upstream_tool=vcflib vcfgeno2haplo
upstream_version=$vcflib_version
fixture=$(realpath --relative-to="$repo_root" "$vcf" 2>/dev/null || printf '%s' "$vcf")
reference=$(realpath --relative-to="$repo_root" "$ref" 2>/dev/null || printf '%s' "$ref")
sample=$sample
window=$window
compared_projection=CHROM,POS,REF,ALT,GT
note=semantic projection only; outputs are not expected to be byte-identical
EOF
