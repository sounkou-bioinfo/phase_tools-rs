#!/usr/bin/env bash
set -euo pipefail

# Default to the repository's explicit public fixture.  For a private/local
# dataset, set VCF, REF, and SAMPLE explicitly in your shell; no private paths
# are embedded in this repository.
repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
c_bin=${C_BIN:-$repo_root/c/phase_mnv}
rs_bin=${RS_BIN:-$repo_root/target/release/phase_mnv_rs}

if [[ -z "${VCF:-}" && -z "${REF:-}" && -z "${SAMPLE:-}" ]]; then
  exec "$repo_root/tests/byte_identical_synthetic.sh" "$rs_bin" "$c_bin"
fi

: "${VCF:?Set VCF=/path/to/input.vcf[.gz|.bcf] for external byte-identity testing}"
: "${REF:?Set REF=/path/to/reference.fa for external byte-identity testing}"
: "${SAMPLE:?Set SAMPLE=<sample-name> for external byte-identity testing}"

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

"$c_bin" -r "$REF" -s "$SAMPLE" -o "$tmp/c.vcf" "$VCF" 2> "$tmp/c.log"
"$rs_bin" -r "$REF" -s "$SAMPLE" -o "$tmp/rs.vcf" "$VCF" 2> "$tmp/rs.log"

cmp "$tmp/c.vcf" "$tmp/rs.vcf"
cmp "$tmp/c.log" "$tmp/rs.log"

echo "phase_mnv_rs byte-identical to C for $VCF"
