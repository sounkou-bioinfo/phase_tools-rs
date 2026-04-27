#!/usr/bin/env bash
set -euo pipefail

case "$(uname -s)" in
  Linux)
    target=${STATIC_TARGET:-$(rustc -vV | sed -n 's/^host: //p')}
    RUSTFLAGS=${STATIC_RUSTFLAGS:-"-C target-feature=+crt-static"} \
      cargo build --release --target "$target"
    ;;
  Darwin)
    echo "macOS does not support fully static executables; building bundled-htslib release." >&2
    if [[ -n "${TARGET:-}" ]]; then
      cargo build --release --target "$TARGET"
    else
      cargo build --release
    fi
    ;;
  *)
    if [[ -n "${TARGET:-}" ]]; then
      cargo build --release --target "$TARGET"
    else
      cargo build --release
    fi
    ;;
esac
