#!/usr/bin/env bash
set -euo pipefail

bin=${1:?usage: $0 <fermi_lite_assemble_binary>}
case "$(uname -m)" in
  x86_64|amd64) ;;
  *) echo "fermi-lite FFI smoke test skipped on non-x86_64 host"; exit 0 ;;
esac

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

"$bin" \
  --min-asm-ovlp 12 \
  --min-count 1 \
  --max-count 1000 \
  --seq ACGTGCTAGCTAGGATCCGATCGATGCTAGCTAGCATC \
  --seq GCTAGGATCCGATCGATGCTAGCTAGCATCGATCGTTA \
  --seq GATCGATGCTAGCTAGCATCGATCGTTAGCTAGGCTAA \
  --seq GCTAGCATCGATCGTTAGCTAGGCTAACCGGTTAAC \
  --seq CGATCGTTAGCTAGGCTAACCGGTTAACCGGTTAGC \
  --seq GCTAGGCTAACCGGTTAACCGGTTAGCATCGATCG \
  > "$tmp/unitigs.fa"

grep -q '^>utg1 ' "$tmp/unitigs.fa"
awk 'BEGIN {seq=0} /^>/ {next} {seq += length($0)} END {exit !(seq >= 40)}' "$tmp/unitigs.fa"

cat > "$tmp/reads.txt" <<'EOF'
>ignored_header
ACGTGCTAGCTAGGATCCGATCGATGCTAGCTAGCATC
GCTAGGATCCGATCGATGCTAGCTAGCATCGATCGTTA
GATCGATGCTAGCTAGCATCGATCGTTAGCTAGGCTAA
GCTAGCATCGATCGTTAGCTAGGCTAACCGGTTAAC
CGATCGTTAGCTAGGCTAACCGGTTAACCGGTTAGC
GCTAGGCTAACCGGTTAACCGGTTAGCATCGATCG
EOF
"$bin" --min-asm-ovlp 12 --min-count 1 --max-count 1000 < "$tmp/reads.txt" > "$tmp/stdin.fa"
grep -q '^>utg1 ' "$tmp/stdin.fa"

echo "fermi-lite FFI smoke test passed"
