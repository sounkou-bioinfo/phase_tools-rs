#!/usr/bin/env bash
set -euo pipefail

bin=${1:?usage: $0 <phase_mnv_binary>}
repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
fixtures="$repo_root/tests/fixtures"
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

cp "$fixtures/ref.fa" "$tmp/ref.fa"
ref="$tmp/ref.fa"

run_body() {
  local input=$1
  local output=$2
  shift 2
  "$bin" -q -r "$ref" -s S1 "$@" "$input" > "$tmp/out.vcf"
  grep -v '^#' "$tmp/out.vcf" > "$output" || true
}

run_body "$fixtures/phased_mnv.vcf" "$tmp/phased.body"
diff -u "$fixtures/phased_mnv.expected.body.vcf" "$tmp/phased.body"

if command -v bcftools >/dev/null 2>&1; then
  bcftools view -Ob -o "$tmp/phased.bcf" "$fixtures/phased_mnv.vcf"
  run_body "$tmp/phased.bcf" "$tmp/phased_bcf.body"
  diff -u "$fixtures/phased_mnv.expected.body.vcf" "$tmp/phased_bcf.body"
fi

run_body "$fixtures/gap.vcf" "$tmp/gap0.body"
diff -u "$fixtures/gap.max0.expected.body.vcf" "$tmp/gap0.body"

run_body "$fixtures/gap.vcf" "$tmp/gap1.body" --max-gap 1
diff -u "$fixtures/gap.max1.expected.body.vcf" "$tmp/gap1.body"

run_body "$fixtures/complex.vcf" "$tmp/complex.body"
diff -u "$fixtures/complex.expected.body.vcf" "$tmp/complex.body"

echo "phase_mnv tests passed"
