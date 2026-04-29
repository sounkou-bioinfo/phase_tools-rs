#!/usr/bin/env bash
set -euo pipefail

# Download Ensembl GRCh37 CDS annotation and build the BED-like codon map used
# by `phase_mnv_rs --mnv-algorithm nirvana-codon`.
#
# The generated files can be large and are intentionally ignored by git.

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
out_dir=${OUT_DIR:-"$repo_root/resources"}
release=${ENSEMBL_GRCH37_RELEASE:-87}
url=${GTF_URL:-"https://ftp.ensembl.org/pub/grch37/release-${release}/gtf/homo_sapiens/Homo_sapiens.GRCh37.${release}.gtf.gz"}
gtf="$out_dir/Homo_sapiens.GRCh37.${release}.gtf.gz"
positions_vcf=${VCF:-}
if [[ -n "$positions_vcf" ]]; then
  default_codon_map="$out_dir/grch37.ensembl${release}.codons.for-vcf.tsv"
else
  default_codon_map="$out_dir/grch37.ensembl${release}.codons.tsv"
fi
codon_map=${CODON_MAP:-"$default_codon_map"}

mkdir -p "$out_dir"

if [[ ! -s "$gtf" ]]; then
  echo "download_grch37_codon_map: downloading $url" >&2
  if command -v curl >/dev/null 2>&1; then
    curl -L --fail --retry 3 -o "$gtf.tmp" "$url"
  elif command -v wget >/dev/null 2>&1; then
    wget -O "$gtf.tmp" "$url"
  else
    echo "error: curl or wget is required" >&2
    exit 1
  fi
  mv "$gtf.tmp" "$gtf"
fi

echo "download_grch37_codon_map: building $codon_map" >&2
if [[ -n "$positions_vcf" ]]; then
  python3 "$repo_root/scripts/gtf_to_codon_map.py" "$gtf" --positions-vcf "$positions_vcf" > "$codon_map.tmp"
else
  python3 "$repo_root/scripts/gtf_to_codon_map.py" "$gtf" > "$codon_map.tmp"
fi
mv "$codon_map.tmp" "$codon_map"

cat <<EOF
codon_map=$codon_map
gtf=$gtf
source_url=$url
positions_vcf=${positions_vcf:-}
EOF
