# fermi-lite vendored source

This directory contains vendored fermi-lite source code from:

- Repository: https://github.com/lh3/fermi-lite
- Upstream commit used for vendoring: `85f159e5a2e15a4c4fb41cdfce91ba0c9b68d2e1`
- License: MIT, see `LICENSE.txt`

The Rust crate builds the library sources through `build.rs` and feeds local
read sequences from Rust. FFI bindings are generated from `fml.h` with bindgen
when libclang is available; `src/fermi_lite_bindings.rs` is a narrow checked-in
fallback for build environments without libclang. Set `PHASE_MNV_REQUIRE_BINDGEN=1`
to fail instead of falling back. The wrapper is currently enabled on x86_64
only because the vendored upstream `ksw.c` path requires SSE2. `bseq.c` is kept
for source attribution/completeness but is not compiled; `src/fermi_lite_shim.c` supplies
the two small sequence helper functions needed by the compiled assembler objects
without linking zlib FASTA/FASTQ input helpers.

If fermi-lite-backed assembly results are used in analyses, cite the FermiKit
paper recommended by upstream fermi-lite:

Li H. FermiKit: assembly-based variant calling for Illumina resequencing data.
Bioinformatics. 2015;31(22):3694-3696. doi:10.1093/bioinformatics/btv440.
