# phase_mnv_rs

Phased VCF/BCF haplotype merger that emits normalized `MNV` or `COMPLEX`
records from variants carried on the same phased haplotype.

The repository contains two implementations:

- `src/main.rs`: Rust implementation, depending on the `rust-htslib` crate.
- `c/phase_mnv.c`: C/htslib reference implementation kept in-tree for
  validation and byte-identical regression tests.

The Rust and C outputs are expected to be byte-identical for the supported
scope. The Rust implementation depends on `rust-htslib` instead of directly on
`hts-sys` so it can use higher-level BAM/BCF APIs for experimental
WhatsHap-inspired read-backed phasing while still dropping to
`rust_htslib::htslib` where exact htslib behavior is needed.

## Citation

This tool's output normalization deliberately implements the left-aligned and
parsimonious variant representation defined by:

> Tan A, Abecasis GR, Kang HM. **Unified representation of genetic variants.**
> *Bioinformatics.* 2015;31(13):2202-2204.
> [doi:10.1093/bioinformatics/btv112](https://doi.org/10.1093/bioinformatics/btv112).

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

## CLI help

The help text below is generated directly from the built binaries when
`README.md` is rendered from `README.Rmd`.

### Rust binary (`phase_mnv_rs`)

```text
usage: phase_mnv -r ref.fa [options] input.vcf|input.bcf

Build minimal merged MNV/complex records from phased variants in one sample.

required:
  -r, --reference FILE   Indexed or indexable FASTA reference

options:
  -s, --sample NAME      Sample to read (default: first sample)
  -o, --output FILE      Output VCF path (default: stdout; plain text)
  -g, --max-gap N        Allow up to N unchanged reference bases between
                        phased variants when building one merged call (default: 0)
      --min-vars N       Minimum source variants per emitted call (default: 2)
      --min-snvs N       Alias for --min-vars
      --unsupported-alleles MODE
                        Selected unsupported allele policy: skip or fail
                        (default: skip)
      --phase-from-bam FILE
                        Experimental Rust read-backed phasing from indexed BAM/CRAM
                        before MNV construction; input GT phase/PS is ignored
      --phase-min-mapq N  Minimum read MAPQ for --phase-from-bam (default: 20)
      --phase-min-baseq N Minimum base quality for --phase-from-bam (default: 13)
      --warn-on-n        Warn when a selected REF/ALT allele contains N
      --no-ref-check     Do not fail when VCF REF differs from FASTA
      --no-header        Suppress VCF header
  -q, --quiet            Suppress summary on stderr
  -h, --help             Show this help

Notes:
  * Only phased diploid GT (e.g. 0|1, 1|0, 1|1, 1|2) is used.
    Unphased, missing, and non-diploid genotypes are skipped.
  * Multi-allelic input sites use the ALT allele selected by each
    haplotype's GT allele index; unselected ALTs are ignored and output
    remains biallelic. Example: GT 1|2 uses ALT1 on haplotype 1 and
    ALT2 on haplotype 2.
  * Symbolic, breakend, spanning-deletion '*', and non-DNA ALT alleles
    are skipped by default and are currently not barriers; use
    --unsupported-alleles fail to reject selected unsupported alleles.
  * FORMAT/PS is honored when present; variants are only merged within the
    same phase set. If PS is absent, the phase separator and proximity
    define the merge block.
  * --phase-from-bam is a Rust-only experimental phaser inspired by
    WhatsHap's read-backed phasing model. It currently phases variants by
    read-supported allele co-occurrence in connected components.
  * With the default --max-gap 0, only adjacent phased variants are
    merged. Pure SNV blocks are TYPE=MNV; blocks containing indels are
    TYPE=COMPLEX.
  * Unless --quiet is set, summary stats go to stderr and include
    input/reference/output (output=stdout for VCF stdout), settings,
    skip counts, unsupported categories, and N counts.
```

### C binary (`phase_mnv`)

```text
usage: phase_mnv -r ref.fa [options] input.vcf|input.bcf

Build minimal merged MNV/complex records from phased variants in one sample.

required:
  -r, --reference FILE   Indexed or indexable FASTA reference

options:
  -s, --sample NAME      Sample to read (default: first sample)
  -o, --output FILE      Output VCF path (default: stdout; plain text)
  -g, --max-gap N        Allow up to N unchanged reference bases between
                        phased variants when building one merged call (default: 0)
      --min-vars N       Minimum source variants per emitted call (default: 2)
      --min-snvs N       Alias for --min-vars
      --unsupported-alleles MODE
                        Selected unsupported allele policy: skip or fail
                        (default: skip)
      --warn-on-n        Warn when a selected REF/ALT allele contains N
      --no-ref-check     Do not fail when VCF REF differs from FASTA
      --no-header        Suppress VCF header
  -q, --quiet            Suppress summary on stderr
  -h, --help             Show this help

Notes:
  * Only phased diploid GT (e.g. 0|1, 1|0, 1|1, 1|2) is used.
    Unphased, missing, and non-diploid genotypes are skipped.
  * Multi-allelic input sites use the ALT allele selected by each
    haplotype's GT allele index; unselected ALTs are ignored and output
    remains biallelic. Example: GT 1|2 uses ALT1 on haplotype 1 and
    ALT2 on haplotype 2.
  * Symbolic, breakend, spanning-deletion '*', and non-DNA ALT alleles
    are skipped by default and are currently not barriers; use
    --unsupported-alleles fail to reject selected unsupported alleles.
  * FORMAT/PS is honored when present; variants are only merged within the
    same phase set. If PS is absent, the phase separator and proximity
    define the merge block.
  * With the default --max-gap 0, only adjacent phased variants are
    merged. Pure SNV blocks are TYPE=MNV; blocks containing indels are
    TYPE=COMPLEX.
  * Unless --quiet is set, summary stats go to stderr and include
    input/reference/output (output=stdout for VCF stdout), settings,
    skip counts, unsupported categories, and N counts.
```

## Examples

These examples are executed when `README.md` is rendered from `README.Rmd`, so
shown output is produced by the current binary instead of copied by hand.

### Adjacent phased SNVs become `TYPE=MNV`

Command:

```bash
target/release/phase_mnv_rs -r tests/fixtures/ref.fa --no-header -q tests/fixtures/phased_mnv.vcf
```

Output:

```text
chr1	1	.	AC	GT	.	PASS	TYPE=MNV;NVAR=2;NSNPS=2;END=2;SOURCE_POS=1,2;HAPS=2;PS=10	GT:PS	0|1:10
chr1	4	.	TA	CG	.	PASS	TYPE=MNV;NVAR=2;NSNPS=2;END=5;SOURCE_POS=4,5;HAPS=1;PS=10	GT:PS	1|0:10
chr1	6	.	CG	AT	.	PASS	TYPE=MNV;NVAR=2;NSNPS=2;END=7;SOURCE_POS=6,7;HAPS=1,2;PS=20	GT:PS	1|1:20
```

### Mixed SNV/indel blocks become `TYPE=COMPLEX`

Command:

```bash
target/release/phase_mnv_rs -r tests/fixtures/ref.fa --no-header -q tests/fixtures/complex.vcf
```

Output:

```text
chr1	2	.	C	TG	.	PASS	TYPE=COMPLEX;NVAR=2;NSNPS=1;END=2;SOURCE_POS=1,2;HAPS=1;PS=30	GT:PS	1|0:30
```

### Rust/C byte-identity smoke test

Command:

```bash
./tests/byte_identical_synthetic.sh target/release/phase_mnv_rs c/phase_mnv
```

Output:

```text
Rust and C outputs are byte-identical on explicit fixture tests/fixtures/byte_identity.vcf
```

### Experimental Rust BAM-backed phasing before MNV construction

The Rust binary can now phase from an indexed BAM/CRAM before constructing MNVs:

```bash
target/release/phase_mnv_rs \
  --reference ref.fa \
  --sample S1 \
  --phase-from-bam reads.bam \
  --max-gap 100 \
  --output mnv.vcf \
  chromosome.vcf.gz
```

This Rust-only mode ignores input `GT` phase and `FORMAT/PS`, extracts allele
co-occurrence from reads with `rust-htslib`, builds connected phase components,
assigns deterministic `PS` values from the first variant in each component, and
then runs the normal MNV/COMPLEX construction. It is intentionally described as
experimental: it is WhatsHap-inspired, but not yet a full clone of WhatsHap's
PedMEC/MEC optimization.

### BAM-backed phasing with external WhatsHap before MNV construction

For the established upstream phaser, or for comparison against the Rust phaser,
use the helper workflow below. It first converts all GT separators from `|` to
`/`, drops `FORMAT/PS` and `FORMAT/PQ` by default, runs `whatshap phase`, then
runs `phase_mnv_rs` on the WhatsHap-phased VCF.

Install the external phasing/indexing tools if needed:

```bash
micromamba create -y -n phase-mnv-whatshap -c conda-forge -c bioconda \
  whatshap bcftools samtools
micromamba run -n phase-mnv-whatshap whatshap --version
```

Run the local BAM-backed phasing workflow with explicit local paths:

```bash
./scripts/phase_from_bam_then_mnv.sh \
  --reference ref.fa \
  --bam reads.bam \
  --vcf chromosome.vcf.gz \
  --sample S1 \
  --max-gap 100 \
  --out-dir local_runs/S1-chromosome
```

The helper writes:

```text
PREFIX.unphased.vcf.gz
PREFIX.whatshap.vcf.gz
PREFIX.phase_mnv.vcf
PREFIX.whatshap.log
PREFIX.phase_mnv.log
```

The unphasing step can also be run directly:

```bash
python3 scripts/unphase_vcf.py input.vcf.gz | bgzip -c > input.unphased.vcf.gz
bcftools index -f input.unphased.vcf.gz
```

### Optional local larger VCF example

The committed README only runs tracked fixtures. To exercise a larger local VCF
while rendering this README, provide explicit paths via environment variables;
the rendered external-example output intentionally omits those paths.

```bash
PHASE_MNV_EXAMPLE_VCF=input.vcf.gz \
PHASE_MNV_EXAMPLE_REF=ref.fa \
PHASE_MNV_EXAMPLE_SAMPLE=S1 \
make readme-external-example
```

Set `PHASE_MNV_EXAMPLE_VCF`, `PHASE_MNV_EXAMPLE_REF`, and
`PHASE_MNV_EXAMPLE_SAMPLE`, then run `make readme-external-example`
to run a local larger VCF through both binaries.

This is intentionally disabled in the committed README so no local/private
paths or data names are embedded.

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
tests/fixtures/multiallelic.vcf
tests/fixtures/multiallelic.expected.body.vcf
tests/fixtures/symbolic.vcf
tests/fixtures/symbolic.max1.expected.body.vcf
tests/fixtures/n_base.vcf
tests/fixtures/n_base.expected.body.vcf
tests/fixtures/read_phase.vcf
tests/fixtures/read_phase.sam
tests/fixtures/read_phase.bam
tests/fixtures/read_phase.bam.bai
tests/fixtures/read_phase.expected.body.vcf
tests/fixtures/byte_identity.vcf
tests/fixtures/ref_mismatch.vcf
tests/fixtures/truncated.vcf.gz
tests/fixtures/vcfgeno2haplo_compare.vcf
```

Run Rust behavior and negative/failure-mode tests:

```bash
make test
```

Run only Rust negative/failure-mode tests:

```bash
make negative-test
```

Run C behavior and negative/failure-mode tests:

```bash
make c-test
```

Run only C negative/failure-mode tests:

```bash
make c-negative-test
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

## README generation

`README.md` is generated from `README.Rmd` so the CLI help and fixture example
outputs stay synchronized with the installed tools:

```bash
make readme
```

This requires R with the `knitr` package. The `readme` target builds both
binaries first, then renders `README.Rmd`. For a local larger VCF example, set
`PHASE_MNV_EXAMPLE_VCF`, `PHASE_MNV_EXAMPLE_REF`, and
`PHASE_MNV_EXAMPLE_SAMPLE`, then run:

```bash
make readme-external-example
```

## Reference implementations

Reference implementations are not vendored in this repository. Clone them when
needed for comparison or compatibility work:

```bash
./scripts/clone_reference_impls.sh
```

### Optional `vcflib vcfgeno2haplo` comparison

`vcfgeno2haplo` is the closest conceptual match for phased genotype-to-haplotype
allele construction, but it is not a byte-identity oracle for this tool. It does
not emit the same `TYPE=MNV`/`TYPE=COMPLEX` schema and it ignores `FORMAT/PS`, so
we compare a narrow semantic projection on a fixture designed for the overlapping
scope: adjacent phased SNVs in one sample.

Run the comparison with micromamba-managed vcflib:

```bash
make compare-vcflib
```

The script uses `vcfgeno2haplo` from `PATH` when available, or creates/runs a
micromamba environment named `phase-mnv-vcflib` with `vcflib=1.0.15` by default.
Override with `VCFLIB_ENV` or `VCFLIB_SPEC` when needed. It compares these
projected columns only:

```text
CHROM POS REF ALT GT
```

Detailed semantics and validation notes are in:

- [`docs/semantics.md`](docs/semantics.md)
- [`docs/validation.md`](docs/validation.md)

The most relevant upstream tools are:

- `vcflib vcfgeno2haplo` — closest conceptual match for converting phased
  genotypes in a window into haplotype alleles.
- `whatshap phase` / `whatshap haplotag` — established read-backed phasing
  reference for future BAM/CRAM work.
- `vt normalize` — normalization reference for left-aligned + parsimonious
  representation.
- `bcftools norm` — useful validator for emitted normalized records.

## CI

GitHub Actions builds and tests both implementations on Linux and macOS:

- Rust release/static where supported
- C binary with bundled static `libhts.a`
- behavior and negative/failure-mode tests for both binaries
- byte-identical Rust-vs-C synthetic fixture test
- Linux semantic-projection comparison against `vcflib=1.0.15` `vcfgeno2haplo`
- binary artifact upload with SHA256 sums

## Notes

- Reads VCF/BCF via htslib/rust-htslib.
- Uses phased diploid GT and `FORMAT/PS`.
- Emits `TYPE=MNV` for pure SNV blocks and `TYPE=COMPLEX` for blocks including indels.
- Normalizes internally with the Tan, Abecasis & Kang 2015 left-aligned + parsimonious rules (doi:10.1093/bioinformatics/btv112).
- Does not require a separate `vt normalize` or `bcftools norm` pass for emitted records.
