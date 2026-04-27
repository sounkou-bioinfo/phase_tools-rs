#!/usr/bin/env bash
set -euo pipefail

rs_bin=${1:?usage: $0 <rust-bin> <c-bin>}
c_bin=${2:?usage: $0 <rust-bin> <c-bin>}
repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
fixtures="$repo_root/tests/fixtures"

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

cp "$fixtures/ref.fa" "$tmp/ref.fa"

"$c_bin" -r "$tmp/ref.fa" -s S1 -o "$tmp/c.vcf" "$fixtures/byte_identity.vcf" 2> "$tmp/c.log"
"$rs_bin" -r "$tmp/ref.fa" -s S1 -o "$tmp/rs.vcf" "$fixtures/byte_identity.vcf" 2> "$tmp/rs.log"

cmp "$tmp/c.vcf" "$tmp/rs.vcf"
cmp "$tmp/c.log" "$tmp/rs.log"

echo "Rust and C outputs are byte-identical on explicit fixture tests/fixtures/byte_identity.vcf"
