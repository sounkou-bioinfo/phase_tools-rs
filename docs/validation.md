# Validation notes

This Rust-only repository validates `phase_tools-rs` binaries with explicit
tracked fixtures. No private/local paths are embedded in tests or documentation.

## Positive behavior fixtures

Run:

```bash
make test
```

Behavior fixtures cover:

- adjacent phased SNVs that become `TYPE=MNV`
- `--max-gap` behavior
- mixed SNV/indel blocks that become `TYPE=COMPLEX`
- selected ALT handling at multi-allelic sites
- symbolic/non-DNA ALT skipping semantics
- `--warn-on-n` warnings while preserving `N` as plain DNA
- native `unphase_vcf` conversion of phased GT separators to unphased GT
  while dropping phase-specific FORMAT tags by default on VCF, stdin VCF, BCF,
  and BGZF VCF output paths
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
- bindgen-backed fermi-lite FFI smoke coverage through `fermi_lite_assemble`,
  including FASTQ/base-quality passthrough
- empirical BAM error-model summary, exact-Q event composition TSV,
  mapchk-like high-nonref site guard, and per-read-position TSV coverage through
  `bam_error_model` on tracked BAM/SAM fixtures, with no MAPQ filter by default
- experimental `phase_adjudicate` pair-level read-evidence coverage on tracked
  VCF/BAM fixtures, including same-phase rows, switched truth/query parity, and
  explicit baseQ filtering
- experimental `bam_contamination` anchor-site contamination probe coverage on
  tracked BAM fixtures, including homozygous-alt reference infiltration,
  optional CHARR-like allele-frequency adjustment, explicit baseQ filtering, and
  unsupported genotype rejection

## Negative/failure-mode fixtures

Run:

```bash
make negative-test
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

1. unphase the input VCF with the native `unphase_vcf` binary;
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
pairs, switch errors, switch rate, blockwise Hamming distance, and blockwise
Hamming rate.

Important limitation: `phase_compare` is not a generic hap.py replacement. It is
for exact-site phasing/block concordance after both paths have been normalized to
the same input records. It does not perform variant representation matching,
ROC/stratification, decompose/atomize, or truth-query callset scoring.

The comparison script accepts thresholds by environment variable:

```bash
MAX_SWITCH_ERRORS=0 MAX_SWITCH_RATE=0 make compare-whatshap-phase
```

For exploratory local runs where non-perfect concordance is expected:

```bash
ALLOW_NONPERFECT=1 KEEP_TMP=1 make compare-whatshap-phase
```

The script sanitizes generated VCF headers before comparison to remove command
lines and local path-bearing records.

## Local/private data policy

Default tests must use tracked fixtures only. Larger validation runs should be
launched with environment overrides, for example:

```bash
WHATSHAP_BIN=whatshap make compare-whatshap-phase
WHATSHAP_ENV=my-whatshap-env make compare-whatshap-phase
REF=ref.fa VCF=input.vcf.gz BAM=reads.bam SAMPLE=S1 ALLOW_NONPERFECT=1 make compare-whatshap-phase
```

Do not commit private paths, sample names, references, BAMs, or generated local
outputs. Use ignored directories such as `local_runs/` or `resources/` for local
experiments.
