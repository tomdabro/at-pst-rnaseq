#!/usr/bin/env bash
# Build a patched STAR for macOS (usage: build_star_macos.sh <output binary>).
#
# The bioconda/upstream macOS STAR maps 0 reads: it aliases its chunk buffers with
# std::stringbuf::pubsetbuf(), which libc++ (macOS) ignores. The patch replaces that
# with an explicit streambuf. Requires Xcode command line tools (clang++, make, curl)
# and libomp (`brew install libomp`).
set -euo pipefail

mkdir -p "$(dirname "$1")"
out="$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
version="2.7.11b"
here="$(cd "$(dirname "$0")" && pwd)"
patch_file="$here/../patches/star-libcxx-streambuf.patch"

omp="$(brew --prefix libomp 2>/dev/null || true)"
if [ -z "$omp" ] || [ ! -d "$omp/include" ]; then
    echo "libomp not found: run 'brew install libomp'" >&2
    exit 1
fi

work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT
curl -sSfL "https://github.com/alexdobin/STAR/archive/refs/tags/$version.tar.gz" | tar xz -C "$work"
cd "$work/STAR-$version"
patch -p1 < "$patch_file"
cd source

# CXXFLAGS_SIMD="": opal's AVX2 path is only used for adapter clipping (--clip3pAdapterSeq),
# which this workflow does not use, and -mavx2 is rejected on arm64.
# LDFLAGS_Mac_static: drop -static-libgcc, which clang rejects.
make -j4 STARforMacStatic CXX=clang++ CXXFLAGS_SIMD="" \
    CXXFLAGS_common="-std=c++11 -Xpreprocessor -fopenmp -I$omp/include "'$(COMPTIMEPLACE) $(GIT_BRANCH_COMMIT_DIFF)' \
    LDFLAGS_Mac_static="-pthread -lz htslib/libhts.a" \
    LDFLAGSextra="-L$omp/lib -lomp" > "$work/build.log" 2>&1 || { tail -30 "$work/build.log" >&2; exit 1; }

cp STAR "$out"
