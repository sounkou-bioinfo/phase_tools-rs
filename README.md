# phase_mnv_rs

Phased VCF/BCF haplotype merger that emits normalized `MNV` or `COMPLEX`
records from variants carried on the same phased haplotype.

The repository contains two implementations:

- `src/main.rs`: Rust implementation, depending on the `rust-htslib` crate.
- `c/phase_mnv.c`: C/htslib reference implementation kept in-tree for
  validation and byte-identical regression tests.

The Rust and C outputs are expected to be byte-identical for the supported
scope. The Rust implementation depends on `rust-htslib` instead of directly on
`hts-sys` so future work can reuse higher-level BAM/BCF APIs for Whatshap-like
read-backed phasing while still dropping to `rust_htslib::htslib` where exact
htslib behavior is needed.

## Citation

This tool's output normalization deliberately implements the left-aligned and
parsimonious variant representation defined by:

> Tan A, Abecasis GR, Kang HM. **Unified representation of genetic variants.**
> *Bioinformatics.* 2015;31(13):2202-2204.
> doi:[10.1093/bioinformatics/btv112](https://doi.org/10.1093/bioinformatics/btv112).

Please cite that paper when relying on phase_mnv_rs normalized `REF`/`ALT`
representation. The emitted VCF header also records the citation in
`##phase_mnv_normalization_citation`.

## Build

```bash
. "$HOME/.cargo/env"
make release
```

The normal Rust release binary is:

```text
target/release/phase_mnv_rs
```

For a fully static Linux Rust binary (static PIE; no `libhts.so`, no glibc
runtime dependency):

```bash
make static-release
host=$(rustc -vV | sed -n 's/^host: //p')
ldd target/$host/release/phase_mnv_rs  # should say statically linked on Linux
```

The Makefile does not hardcode a Rust target. Override when needed:

```bash
make release TARGET=aarch64-apple-darwin
make static-release STATIC_TARGET=x86_64-unknown-linux-gnu
```

Install to `~/.local/bin`:

```bash
make install          # normal release
make install-static   # static release on Linux, bundled-htslib release elsewhere
```

## C implementation

Build against system htslib:

```bash
make -C c
```

Build a bundled/static-htslib C binary:

```bash
make c-static
```

On Linux this script attempts a fully static executable. On macOS, fully static
executables are not supported by the platform; the script links `libhts.a`
statically while system libraries remain dynamic.

## Test

All CI tests use explicit files under `tests/fixtures/`:

```text
tests/fixtures/ref.fa
tests/fixtures/phased_mnv.vcf
tests/fixtures/phased_mnv.expected.body.vcf
tests/fixtures/gap.vcf
tests/fixtures/gap.max0.expected.body.vcf
tests/fixtures/gap.max1.expected.body.vcf
tests/fixtures/complex.vcf
tests/fixtures/complex.expected.body.vcf
tests/fixtures/byte_identity.vcf
```

Run Rust behavior tests:

```bash
make test
```

Run C behavior tests:

```bash
make c-test
```

Compare Rust and C byte-for-byte on the explicit byte-identity fixture:

```bash
./tests/byte_identical_synthetic.sh target/release/phase_mnv_rs c/phase_mnv
```

Compare byte-for-byte against the C tool on the explicit public fixture:

```bash
make byte-test
```

For private/local datasets, provide paths explicitly from your shell; no private
paths are embedded in this repository:

```bash
VCF=input.vcf.gz REF=ref.fa SAMPLE=S1 make byte-test
```

## CI

GitHub Actions builds and tests both implementations on Linux and macOS:

- Rust release/static where supported
- C binary with bundled static `libhts.a`
- behavior tests for both binaries
- byte-identical Rust-vs-C synthetic fixture test
- binary artifact upload with SHA256 sums

## Notes

- Reads VCF/BCF via htslib/rust-htslib.
- Uses phased diploid GT and `FORMAT/PS`.
- Emits `TYPE=MNV` for pure SNV blocks and `TYPE=COMPLEX` for blocks including indels.
- Normalizes internally with the Tan, Abecasis & Kang 2015 left-aligned + parsimonious rules (doi:10.1093/bioinformatics/btv112).
- Does not require a separate `vt normalize` or `bcftools norm` pass for emitted records.
