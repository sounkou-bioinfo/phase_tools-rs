#!/usr/bin/env bash
set -euo pipefail

root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
htslib_version=${HTSLIB_VERSION:-1.19.1}
build_dir=${BUILD_DIR:-"$root/.build/c-static"}
prefix=${HTSLIB_PREFIX:-"$build_dir/htslib-install"}
jobs=${JOBS:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 2)}

mkdir -p "$build_dir"

archive="$build_dir/htslib-${htslib_version}.tar.bz2"
src="$build_dir/htslib-${htslib_version}"
url="https://github.com/samtools/htslib/releases/download/${htslib_version}/htslib-${htslib_version}.tar.bz2"

if [[ ! -f "$prefix/lib/libhts.a" ]]; then
  if [[ ! -f "$archive" ]]; then
    curl -L --retry 3 -o "$archive" "$url"
  fi
  rm -rf "$src"
  tar -xjf "$archive" -C "$build_dir"
  (
    cd "$src"
    ./configure \
      --prefix="$prefix" \
      --disable-bz2 \
      --disable-lzma \
      --disable-libcurl \
      --without-libdeflate
    make -j"$jobs" lib-static
    make install
  )
fi

cc_bin=${CC:-cc}
cflags=${CFLAGS:-"-O3 -Wall -Wextra -std=c11"}
cppflags=${CPPFLAGS:-"-I$prefix/include"}
libs="$prefix/lib/libhts.a -lz -lm -lpthread"
extra_ldflags=""

case "$(uname -s)" in
  Linux)
    # Fully static Linux executable when static zlib/libpthread/libm are present.
    extra_ldflags="-static"
    ;;
  Darwin)
    # macOS does not support fully static executables. This still links libhts.a
    # statically; system libraries remain dynamically linked as required by macOS.
    ;;
  *)
    ;;
esac

# Intentionally use word splitting for user-overridable CFLAGS/CPPFLAGS/LDFLAGS.
# Avoid an empty bash array here: GitHub's macOS bash 3.2 treats "${empty[@]}"
# as an unbound variable under `set -u`.
# shellcheck disable=SC2086
"$cc_bin" $cflags $cppflags "$root/c/phase_mnv.c" $extra_ldflags $libs -o "$root/c/phase_mnv"

"$root/c/phase_mnv" --help >/dev/null
file "$root/c/phase_mnv"
