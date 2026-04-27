# Validation notes

This repository validates `phase_mnv_rs` with explicit tracked fixtures and an
in-tree C/htslib reference implementation. No private/local paths are embedded
in tests or documentation.

## Rust/C byte identity

The strongest compatibility target is byte-for-byte identity between the Rust
implementation and the C implementation for the supported scope.

Default command:

```bash
make byte-test
```

This runs both binaries on:

```text
tests/fixtures/byte_identity.vcf
tests/fixtures/ref.fa
```

and compares both VCF output and stderr logs with `cmp`.

## Positive behavior fixtures

Run:

```bash
make test
make c-test
```

Behavior fixtures cover:

- adjacent phased SNVs that become `TYPE=MNV`
- `--max-gap` behavior
- mixed SNV/indel blocks that become `TYPE=COMPLEX`
- Rust/C byte identity for supported synthetic cases

## Negative/failure-mode fixtures

Run:

```bash
make negative-test
make c-negative-test
```

Negative fixtures and generated checks cover:

- missing required reference option
- missing input VCF/BCF
- missing FASTA reference
- unknown sample
- invalid negative `--max-gap`
- REF/FASTA mismatch
- truncated gzipped VCF input

The truncated-input fixture is tracked explicitly:

```text
tests/fixtures/truncated.vcf.gz
```

## `vcflib vcfgeno2haplo` comparison

`vcflib vcfgeno2haplo` is the closest conceptual upstream tool for converting
phased genotypes within a window into haplotype alleles. It is **not** a
byte-identical oracle for this project because:

- it emits haplotype-allele VCF records, not the `TYPE=MNV` / `TYPE=COMPLEX`
  schema used by `phase_mnv_rs`
- it clusters by a window rule and does not honor `FORMAT/PS` like this tool
- it may pass through non-cluster input records, while this tool emits only
  merged haplotype records
- normalization responsibilities differ

The comparison is therefore intentionally narrow: adjacent phased SNVs in one
sample, projected to fields that both tools can represent directly.

Run:

```bash
make compare-vcflib
```

Default upstream package:

```text
vcflib=1.0.15
```

Default fixture:

```text
tests/fixtures/vcfgeno2haplo_compare.vcf
```

Projection compared:

```text
CHROM POS REF ALT GT
```

The script uses `vcfgeno2haplo` from `PATH` when available. Otherwise it uses
micromamba to create/run an environment named `phase-mnv-vcflib`.

Relevant overrides:

```bash
VCFLIB_ENV=my-env make compare-vcflib
VCFLIB_SPEC='vcflib=1.0.15' make compare-vcflib
VCFGENO2HAPLO_BIN=/path/to/vcfgeno2haplo make compare-vcflib
```

## Normalization references

The emitted `REF`/`ALT` records are normalized internally according to the
left-aligned and parsimonious representation described by:

> Tan A, Abecasis GR, Kang HM. Unified representation of genetic variants.
> Bioinformatics. 2015;31(13):2202-2204. doi:10.1093/bioinformatics/btv112.

`vt normalize` and `bcftools norm` remain useful external validators, but they
are not phase-aware haplotype merging tools and are not required as a post-pass
for supported outputs.
