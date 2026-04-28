#!/usr/bin/env bash
set -euo pipefail

# Clone upstream/reference implementations used for design and validation.
# These repositories are intentionally NOT vendored or committed into this repo.

root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
dest=${REFERENCE_IMPLS_DIR:-"$root/reference_impls"}
mkdir -p "$dest"

clone_or_update() {
  local url=$1
  local name=$2
  if [[ -d "$dest/$name/.git" ]]; then
    git -C "$dest/$name" fetch --tags --prune
  else
    git clone "$url" "$dest/$name"
  fi
  printf '%-10s %s\n' "$name" "$(git -C "$dest/$name" rev-parse HEAD)"
}

clone_or_update https://github.com/vcflib/vcflib.git vcflib
clone_or_update https://github.com/whatshap/whatshap.git whatshap
clone_or_update https://github.com/atks/vt.git vt
clone_or_update https://github.com/samtools/bcftools.git bcftools
clone_or_update https://github.com/Illumina/Nirvana.git Nirvana
