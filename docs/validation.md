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
- selected ALT handling at multi-allelic sites
- symbolic/non-DNA ALT skipping semantics
- `--warn-on-n` warnings while preserving `N` as plain DNA
- `scripts/unphase_vcf.py` conversion of phased GT separators to unphased GT
  while dropping phase-specific FORMAT tags by default
- experimental Rust `--phase-from-bam` read-backed phasing on a tiny tracked
  BAM/BAI fixture before MNV construction, using the default exact
  single-sample MEC dynamic-programming algorithm
- Rust output format inference for plain VCF, BGZF-compressed VCF, and BCF,
  including `--threads` plumbing for compressed input/output checks when
  `bcftools` is available
- `bcftools norm -f ... -c x` checks that emitted MNV/COMPLEX records are not
  further realigned or mismatch-removed on tracked fixtures
- Rust `--emit all-sites` header preservation: the original VCF header is kept
  and `phase_mnv` metadata is appended while BAM-backed `GT:PS` updates are
  applied to all input records in the tiny tracked fixture
- Rust `--mnv-algorithm nirvana-codon` same-codon SNV recomposition on a tracked
  BED-like codon-map fixture
- native `phase_compare` switch-error, phase-match, and blockwise-Hamming stats
  on a tiny tracked truth/query fixture
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
- `--unsupported-alleles fail` on selected unsupported ALT alleles
- truncated gzipped VCF input

The truncated-input fixture is tracked explicitly:

```text
tests/fixtures/truncated.vcf.gz
```

## WhatsHap + native `phase_compare`

The external conformance comparison is WhatsHap-based and no longer uses
hap.py. It uses the in-repository `phase_compare` binary, which is a narrow,
fast phase-concordance comparator.

Run:

```bash
make compare-whatshap-phase
```

The script compares two paths on the tracked tiny BAM/VCF fixture by default:

1. unphase the input VCF with `scripts/unphase_vcf.py`;
2. run external `whatshap phase` on the unphased VCF and BAM to create the truth
   all-sites phased VCF;
3. run Rust `phase_mnv_rs --emit all-sites --phase-from-bam` directly on the
   input VCF and BAM to create the query all-sites phased VCF;
4. run `phase_compare` on truth/query all-sites VCFs.

Default fixture inputs:

```text
tests/fixtures/read_phase.vcf
tests/fixtures/read_phase.bam
tests/fixtures/ref.fa
```

`phase_compare` reports exact shared variant records, common heterozygous sites,
phased sites with PS in both files, intersection PS blocks, assessed adjacent
pairs, phase-match pairs, switch errors, switch rate, and blockwise Hamming
rate. The comparison is a phase/PS-block metric, not a generic variant-call
matching metric.

The generated VCFs are passed through `scripts/sanitize_vcf_headers.py`, which
removes `##bcftools_*` headers, producer command-line headers, and path-bearing
`phase_mnv`/reference header records before comparison. This keeps local
filesystem paths out of comparison VCFs.

The comparison script discovers `whatshap` from `PATH` first, then falls back to
a micromamba environment named `phase-mnv-whatshap`.

Relevant overrides:

```bash
WHATSHAP_BIN=whatshap make compare-whatshap-phase
WHATSHAP_ENV=my-whatshap-env make compare-whatshap-phase
REF=ref.fa VCF=input.vcf.gz BAM=reads.bam SAMPLE=S1 ALLOW_NONPERFECT=1 make compare-whatshap-phase
```

## Nirvana recomposition benchmark target

Illumina Nirvana's MNV recomposition is a better benchmark than generic
haplotype-comparison output when the target is codon/transcript-aware SNV
recomposition. The reference clone helper includes Nirvana, but the repository
is not vendored:

```bash
./scripts/clone_reference_impls.sh
```

The first implemented benchmark slice is `--mnv-algorithm nirvana-codon` (see
`docs/nirvana_benchmark.md`), which uses a small BED-like codon map and emits
SNV-only MNVs where two or more phased SNV observations share a transcript/codon
key. For local GRCh37 experiments, build an ignored Ensembl-derived codon map
with:

```bash
./scripts/download_grch37_codon_map.sh
```

For large real VCFs, prefer a map restricted to SNV positions in that VCF:

```bash
VCF=input.vcf.gz ./scripts/download_grch37_codon_map.sh
```

The broader intended benchmark scope remains Nirvana-like phase-set and
homozygous-variant semantics, adjacent-codon aggregation,
unsupported-overlap barriers, sample-specific multi-sample recomposition, and
exact linkage/quality/filter behavior. Indel/complex recomposition remains a
separate policy decision.

## Normalization references

The emitted `REF`/`ALT` records are normalized internally according to the
left-aligned and parsimonious representation described by:

> Tan A, Abecasis GR, Kang HM. Unified representation of genetic variants.
> Bioinformatics. 2015;31(13):2202-2204. doi:10.1093/bioinformatics/btv112.

`vt normalize` and `bcftools norm` remain useful external validators, but they
are not phase-aware haplotype merging tools and are not required as a post-pass
for supported outputs.
