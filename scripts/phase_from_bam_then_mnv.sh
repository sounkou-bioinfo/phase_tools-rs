#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
usage: phase_from_bam_then_mnv.sh -r ref.fa -b reads.bam -v variants.vcf.gz [options]

Unphase an input VCF, phase it from BAM/CRAM reads with WhatsHap, then run
phase_mnv_rs on the WhatsHap-phased VCF.

required:
  -r, --reference FILE       Reference FASTA used by both WhatsHap and phase_mnv_rs
  -b, --bam FILE             Indexed BAM/CRAM with reads for the selected sample
  -v, --vcf FILE             Input VCF/VCF.GZ/BCF to first make unphased

options:
  -s, --sample NAME          Sample name (default: first VCF sample via bcftools)
  -O, --out-dir DIR          Output directory (default: phase_mnv_from_bam)
  -p, --prefix NAME          Output file prefix (default: sample name)
  -g, --max-gap N            phase_mnv_rs --max-gap value (default: 0)
      --phase-mnv-bin FILE   phase_mnv_rs binary (default: target/release/phase_mnv_rs)
      --unphase-bin FILE     unphase_vcf binary (default: target/release/unphase_vcf)
      --whatshap FILE        whatshap executable (default: whatshap from PATH)
      --whatshap-extra ARGS  Extra shell words appended to `whatshap phase`
      --keep-phase-tags      Keep FORMAT/PS and FORMAT/PQ in the unphased VCF
  -h, --help                 Show this help

outputs in OUT_DIR:
  PREFIX.unphased.vcf.gz       VCF after replacing GT '|' with '/' and dropping PS/PQ
  PREFIX.whatshap.vcf.gz       WhatsHap-phased VCF
  PREFIX.phase_mnv.vcf         phase_mnv_rs output VCF
  PREFIX.phase_mnv.log         phase_mnv_rs stderr summary/statistics
  PREFIX.whatshap.log          WhatsHap stderr/stdout log

No input paths are baked into the repository; provide local BAM/VCF/reference
paths explicitly when running this helper.
USAGE
}

die() {
  printf 'error: %s\n' "$*" >&2
  exit 1
}

script_dir=$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_dir=$(CDPATH= cd -- "$script_dir/.." && pwd)

ref=""
bam=""
vcf=""
sample=""
out_dir="phase_mnv_from_bam"
prefix=""
max_gap="0"
phase_mnv_bin="$repo_dir/target/release/phase_mnv_rs"
unphase_bin="$repo_dir/target/release/unphase_vcf"
whatshap_bin="whatshap"
whatshap_extra=""
keep_phase_tags=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    -r|--reference)
      [[ $# -ge 2 ]] || die "$1 requires an argument"; ref=$2; shift 2 ;;
    -b|--bam|--cram)
      [[ $# -ge 2 ]] || die "$1 requires an argument"; bam=$2; shift 2 ;;
    -v|--vcf)
      [[ $# -ge 2 ]] || die "$1 requires an argument"; vcf=$2; shift 2 ;;
    -s|--sample)
      [[ $# -ge 2 ]] || die "$1 requires an argument"; sample=$2; shift 2 ;;
    -O|--out-dir)
      [[ $# -ge 2 ]] || die "$1 requires an argument"; out_dir=$2; shift 2 ;;
    -p|--prefix)
      [[ $# -ge 2 ]] || die "$1 requires an argument"; prefix=$2; shift 2 ;;
    -g|--max-gap)
      [[ $# -ge 2 ]] || die "$1 requires an argument"; max_gap=$2; shift 2 ;;
    --phase-mnv-bin)
      [[ $# -ge 2 ]] || die "$1 requires an argument"; phase_mnv_bin=$2; shift 2 ;;
    --unphase-bin)
      [[ $# -ge 2 ]] || die "$1 requires an argument"; unphase_bin=$2; shift 2 ;;
    --whatshap)
      [[ $# -ge 2 ]] || die "$1 requires an argument"; whatshap_bin=$2; shift 2 ;;
    --whatshap-extra)
      [[ $# -ge 2 ]] || die "$1 requires an argument"; whatshap_extra=$2; shift 2 ;;
    --keep-phase-tags)
      keep_phase_tags=1; shift ;;
    -h|--help)
      usage; exit 0 ;;
    *)
      usage >&2; die "unknown option: $1" ;;
  esac
done

[[ -n "$ref" ]] || { usage >&2; die "--reference is required"; }
[[ -n "$bam" ]] || { usage >&2; die "--bam is required"; }
[[ -n "$vcf" ]] || { usage >&2; die "--vcf is required"; }
[[ -f "$ref" ]] || die "reference FASTA does not exist: $ref"
[[ -f "$bam" ]] || die "BAM/CRAM does not exist: $bam"
[[ -f "$vcf" ]] || die "VCF/BCF does not exist: $vcf"
[[ "$max_gap" =~ ^[0-9]+$ ]] || die "--max-gap must be a non-negative integer"

command -v bcftools >/dev/null 2>&1 || die "bcftools is required to sample-detect/index outputs"
command -v "$whatshap_bin" >/dev/null 2>&1 || die "whatshap not found: $whatshap_bin"

if [[ ! -x "$phase_mnv_bin" ]]; then
  if [[ "$phase_mnv_bin" == "$repo_dir/target/release/phase_mnv_rs" ]]; then
    (cd "$repo_dir" && cargo build --release)
  fi
fi
[[ -x "$phase_mnv_bin" ]] || die "phase_mnv_rs binary is not executable: $phase_mnv_bin"

if [[ ! -x "$unphase_bin" ]]; then
  if [[ "$unphase_bin" == "$repo_dir/target/release/unphase_vcf" ]]; then
    (cd "$repo_dir" && cargo build --release --bin unphase_vcf)
  fi
fi
[[ -x "$unphase_bin" ]] || die "unphase_vcf binary is not executable: $unphase_bin"

if [[ -z "$sample" ]]; then
  sample=$(bcftools query -l "$vcf" | head -n 1)
  [[ -n "$sample" ]] || die "could not determine sample from VCF; pass --sample"
fi

if [[ -z "$prefix" ]]; then
  prefix=$(printf '%s' "$sample" | sed 's/[^[:alnum:]._-]/_/g')
fi

mkdir -p "$out_dir"
unphased_vcf="$out_dir/$prefix.unphased.vcf.gz"
whatshap_vcf="$out_dir/$prefix.whatshap.vcf.gz"
whatshap_plain="$out_dir/$prefix.whatshap.tmp.vcf"
mnv_vcf="$out_dir/$prefix.phase_mnv.vcf"
whatshap_log="$out_dir/$prefix.whatshap.log"
mnv_log="$out_dir/$prefix.phase_mnv.log"

unphase_args=("$unphase_bin" "$vcf")
if [[ "$keep_phase_tags" -eq 1 ]]; then
  unphase_args+=(--keep-phase-tags)
fi

printf 'phase_from_bam_then_mnv: writing %s\n' "$unphased_vcf" >&2
if command -v bgzip >/dev/null 2>&1; then
  "${unphase_args[@]}" | bgzip -c > "$unphased_vcf"
else
  "${unphase_args[@]}" | bcftools view -Oz -o "$unphased_vcf" -
fi
bcftools index -f "$unphased_vcf"

whatshap_args=(phase --reference "$ref" --sample "$sample" -o "$whatshap_plain")
if [[ -n "$whatshap_extra" ]]; then
  # Intentional shell splitting for advanced local WhatsHap options.
  # shellcheck disable=SC2206
  extra_words=($whatshap_extra)
  whatshap_args+=("${extra_words[@]}")
fi
whatshap_args+=("$unphased_vcf" "$bam")

printf 'phase_from_bam_then_mnv: running WhatsHap -> %s\n' "$whatshap_plain" >&2
"$whatshap_bin" "${whatshap_args[@]}" > "$whatshap_log" 2>&1
printf 'phase_from_bam_then_mnv: compressing WhatsHap VCF -> %s\n' "$whatshap_vcf" >&2
if command -v bgzip >/dev/null 2>&1; then
  bgzip -c "$whatshap_plain" > "$whatshap_vcf"
else
  bcftools view -Oz -o "$whatshap_vcf" "$whatshap_plain"
fi
rm -f "$whatshap_plain"
bcftools index -f "$whatshap_vcf"

printf 'phase_from_bam_then_mnv: running phase_mnv_rs -> %s\n' "$mnv_vcf" >&2
"$phase_mnv_bin" \
  -r "$ref" \
  -s "$sample" \
  --max-gap "$max_gap" \
  -o "$mnv_vcf" \
  "$whatshap_vcf" > "$mnv_log" 2>&1

cat <<EOF
unphased_vcf=$unphased_vcf
whatshap_vcf=$whatshap_vcf
phase_mnv_vcf=$mnv_vcf
whatshap_log=$whatshap_log
phase_mnv_log=$mnv_log
EOF
