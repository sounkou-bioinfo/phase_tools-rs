# phase_mnv_rs

`phase_mnv_rs` is a phased VCF/BCF haplotype merger. It takes variants carried
on the same phased haplotype/phase set and emits normalized `TYPE=MNV` or
`TYPE=COMPLEX` records.

## Scope at a glance

| Area | Status |
| --- | --- |
| Default output | Derived MNV/COMPLEX records only (`--emit mnv`). |
| Canonical implementation | Rust (`src/main.rs`), built on `rust-htslib`. |
| Native dependencies | Vendored C libraries only when useful, exposed through generated/bindgen-backed FFI wrappers. |
| Extra Rust modes | BCF/BGZF output, all-sites output, experimental BAM/CRAM phasing, codon-aware recomposition, native phase comparison, native VCF/BCF unphasing, empirical BAM error-model summaries, and experimental fermi-lite local assembly bindings. |

The Rust binary uses `rust-htslib` for VCF/BCF/BAM/CRAM APIs and only drops to
`rust_htslib::htslib` where exact htslib behavior is needed. Rust-only features
are documented as experimental unless explicitly validated against an upstream
reference.

## Quick start

```bash
. "$HOME/.cargo/env"
make release

target/release/phase_mnv_rs \
  --reference ref.fa \
  --sample S1 \
  --output mnv.vcf \
  input.vcf.gz
```

Build or install variants:

```bash
make install          # install Rust CLIs to ~/.local/bin
make static-release   # static Linux Rust binary when supported
make install-static   # install static Linux binary when supported
# cargo also builds phase_compare, unphase_vcf, fermi_lite_assemble, and bam_error_model
```

Override the Rust target when needed:

```bash
make release TARGET=aarch64-apple-darwin
make static-release STATIC_TARGET=x86_64-unknown-linux-gnu
```

## Modes and output

| Option | Meaning |
| --- | --- |
| `--emit mnv` | Default. Emit only derived `TYPE=MNV` / `TYPE=COMPLEX` records. |
| `--emit all-sites` | Preserve all input records/header and update `GT:PS` where phasing is available. Does not append MNV/COMPLEX records yet. |
| `--mnv-algorithm proximity` | Default window/phase-set recomposition. |
| `--mnv-algorithm nirvana-codon --codon-map codons.tsv` | Initial SNV-only same-codon recomposition slice; not full Nirvana parity. |
| `--phase-from-bam reads.bam` | Experimental indexed BAM/CRAM-backed phasing before output. |

The codon map is BED-like: `CHROM START0 END0 TRANSCRIPT CODON_ID`. In
`nirvana-codon` mode, phased SNV observations are recomposed only when at least
two observations on the same haplotype/phase set share a transcript/codon key.

Output format is inferred from `-o/--output`:

```text
-o out.vcf      plain VCF
-o out.vcf.gz   BGZF-compressed VCF
-o out.vcf.bgz  BGZF-compressed VCF
-o out.bcf      BCF
```

Use `--threads N` or `-@ N` for htslib/BGZF worker threads on compressed input
and output. Stdout remains plain VCF.

## Citation

Please cite the upstream methods that match the parts of `phase_mnv_rs` you use.

For normalized `REF`/`ALT` representation:

> Tan A, Abecasis GR, Kang HM. **Unified representation of genetic variants.**
> *Bioinformatics.* 2015;31(13):2202-2204.
> [doi:10.1093/bioinformatics/btv112](https://doi.org/10.1093/bioinformatics/btv112).

For workflows that run external WhatsHap via
`scripts/phase_from_bam_then_mnv.sh`:

> Martin M, Patterson M, Garg S, Fischer SO, Pisanti N, Klau GW, Schoenhuth A,
> Marschall T. **WhatsHap: fast and accurate read-based phasing.**
> *bioRxiv.* 085050.
> [doi:10.1101/085050](https://doi.org/10.1101/085050).

For the read-backed weighted haplotype assembly/MEC ideas behind the
experimental Rust `--phase-from-bam --phase-algorithm mec` mode:

> Patterson M, Marschall T, Pisanti N, van Iersel L, Stougie L, Klau GW,
> Schönhuth A. **WhatsHap: Weighted Haplotype Assembly for Future-Generation
> Sequencing Reads.** *Journal of Computational Biology.* 2015;22(6):498-509.
> [doi:10.1089/cmb.2014.0157](https://doi.org/10.1089/cmb.2014.0157).

The Rust BAM phaser is WhatsHap-inspired, not a full WhatsHap/PedMEC clone. Use
WhatsHap itself when you need established WhatsHap behavior and cite the
appropriate WhatsHap publication for that workflow.

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
  -o, --output FILE      Output path (default: stdout). Format is inferred:
                        .vcf = plain VCF, .vcf.gz/.vcf.bgz = BGZF VCF,
                        .bcf = BCF; stdout defaults to plain VCF
  -@, --threads N        Extra htslib/BGZF threads for decompression and
                        compressed output (default: 1)
      --emit MODE        Output mode: mnv (default) or all-sites. all-sites
                        is Rust-only, preserves input records/header, and
                        updates GT/PS when used with --phase-from-bam
  -g, --max-gap N        Allow up to N unchanged reference bases between
                        phased variants when building one merged call (default: 0)
      --mnv-algorithm MODE
                        MNV construction: proximity (default) or
                        nirvana-codon (SNV-only same-codon seed mode)
      --codon-map FILE   BED-like codon map for --mnv-algorithm nirvana-codon:
                        CHROM START0 END0 TRANSCRIPT CODON_ID [ignored...]
      --min-vars N       Minimum source variants per emitted call (default: 2)
      --min-snvs N       Alias for --min-vars
      --unsupported-alleles MODE
                        Selected unsupported allele policy: skip or fail
                        (default: skip)
      --phase-from-bam FILE
                        Experimental Rust read-backed phasing from indexed BAM/CRAM
                        before MNV construction; input GT phase/PS is ignored
      --phase-algorithm MODE
                        BAM phasing algorithm: mec or greedy (default: mec)
      --phase-max-coverage N
                        Maximum selected read coverage per variant for MEC phasing
                        (default: 15; WhatsHap-style downsampling guard)
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
    WhatsHap's read-backed phasing model. The default mec algorithm solves
    exact diploid single-sample MEC per connected component after deterministic
    read selection; greedy keeps the earlier pairwise parity heuristic.
  * With the default --max-gap 0 and --mnv-algorithm proximity, only
    adjacent phased variants are merged. Pure SNV blocks are TYPE=MNV;
    blocks containing indels are TYPE=COMPLEX. nirvana-codon mode only
    recomposes phased SNVs sharing a codon key from --codon-map.
  * Output format is inferred from -o/--output. BCF output always includes
    a VCF/BCF header even if --no-header is set.
  * --emit all-sites keeps the original VCF/BCF header via htslib and
    appends phase_mnv metadata instead of replacing it.
  * Unless --quiet is set, summary stats go to stderr and include
    input/reference/output (output=stdout for VCF stdout), settings,
    skip counts, unsupported categories, and N counts.
```

### Phase comparison binary (`phase_compare`)

```text
usage: phase_compare [options] truth.vcf|bcf query.vcf|bcf

Fast phase-concordance comparison for two VCF/BCF files. The tool compares
exact shared variant records, diploid heterozygous GT phase, PS block
intersections, pairwise phase matches, and switch errors. It does not perform
generic haplotype variant-call matching like hap.py.

options:
-s, --sample NAME          Sample name used in both files (default: first truth sample)
--truth-sample NAME    Sample name in truth file
--query-sample NAME    Sample name in query file
--ignore-sample-name   Use the first sample in the query if the truth sample name is absent
--only-snvs            Restrict to heterozygous selected SNV genotypes
-@, --threads N            htslib reader threads (default: 1)
-o, --report-prefix PREFIX Write PREFIX.summary.tsv, like hap.py's -o prefix
--summary-tsv FILE     Write summary TSV to FILE as well as stdout
--switch-bed FILE      Write switch-error intervals as BED
--switch-error-bed FILE Alias for --switch-bed, compatible with whatshap compare
--pair-tsv FILE        Write assessed adjacent-pair decisions
--tsv-pairwise FILE    Alias for --pair-tsv, compatible with whatshap compare
-r, --reference FILE       Accepted for hap.py-style compatibility; ignored
--engine NAME          Accepted for hap.py-style compatibility; ignored
--no-roc               Accepted for hap.py-style compatibility; ignored
--no-decompose         Accepted for hap.py-style compatibility; ignored
--names NAMES          Accepted for whatshap-compare compatibility; ignored
--tsv-multiway FILE    Accepted for whatshap-compare compatibility; ignored
-h, --help                 Show this help

Output is a TSV summary with per-contig rows and a final TOTAL row.
```

### Native unphasing helper (`unphase_vcf`)

```text
usage: unphase_vcf [options] input.vcf|input.vcf.gz|input.bcf|-

Write an unphased VCF stream from VCF/VCF.GZ/BCF input. GT separators are
converted from phased to unphased, and FORMAT/PS plus FORMAT/PQ are removed
by default. Other records, alleles, INFO fields, filters, and non-phase FORMAT
values are preserved through htslib/rust-htslib.

options:
-o, --output FILE       Output VCF path; .gz/.bgz writes BGZF (default: stdout)
--keep-phase-tags   Keep FORMAT/PS and FORMAT/PQ instead of removing them
-@, --threads N         htslib threads for compressed input/output (default: 1)
-h, --help              Show this help
```

### fermi-lite smoke/utility binary (`fermi_lite_assemble`)

```text
usage: fermi_lite_assemble [options] [--seq SEQ ...]

Small fermi-lite FFI smoke/utility binary. With --seq, assembles the supplied
sequences. Without --seq, reads one plain sequence per non-empty stdin line,
ignoring FASTA-style header lines. With --fastq, reads FASTQ from stdin and
passes base qualities to fermi-lite's error-correction path when --ec-k >= 0.
This is intended for local adjudication experiments, not as a full fermi-lite
CLI replacement.

options:
--seq SEQ              Add one input read/sequence
--fastq                Read FASTQ records from stdin instead of plain lines
-@, --threads N            fermi-lite threads (default: 1)
--min-asm-ovlp N       minimum assembly overlap (default: 21)
--min-count N          minimum k-mer count threshold (default: 1)
--max-count N          maximum k-mer count threshold (default: 1000)
--ec-k N               error-correction k; negative disables EC (default: -1)
-h, --help                 Show this help
```

### BAM empirical error model helper (`bam_error_model`)

```text
usage: bam_error_model --reference ref.fa [options] reads.bam|reads.cram

Learn a simple empirical sequencing-error table from aligned reads by comparing
BAM/CRAM bases to a FASTA reference. No MAPQ filter is applied by default; MAPQ
is summarized as a covariate. Known variant sites are not masked yet, so real
biological variants contribute to the mismatch rate unless callers restrict the
regions accordingly.

options:
-r, --reference FILE       FASTA reference with .fai
--region REG          Restrict to region CHR:START-END (1-based, repeatable)
--max-reads N         Stop after N usable reads
--min-mapq N          Optional MAPQ cutoff (default: 0; no cutoff)
--position-tsv FILE   Write per-read-position empirical error TSV
--high-quality-threshold N  BaseQ threshold for high/low position groups (default: 20)
--include-duplicates  Include duplicate reads
--include-secondary   Include secondary alignments
--include-supplementary Include supplementary alignments
-h, --help                Show this help
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
observations from reads with `rust-htslib`, performs deterministic read selection
(default `--phase-max-coverage 15`), builds connected phase components, and by
default solves the exact weighted diploid single-sample MEC objective per
component (`--phase-algorithm mec`). The earlier pairwise parity heuristic is
still available as `--phase-algorithm greedy`. It is intentionally described as
experimental: it is WhatsHap-inspired, but not yet a full clone of WhatsHap's
PedMEC implementation or all WhatsHap options.

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
target/release/unphase_vcf input.vcf.gz | bgzip -c > input.unphased.vcf.gz
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
to run a local larger VCF through the Rust binary.

This is intentionally disabled in the committed README so no local/private
paths or data names are embedded.

## Test

CI and default tests use only tracked fixtures under `tests/fixtures/`.

```bash
make test                  # Rust behavior, unphase_vcf, output-format, bcftools norm, BAM phasing, phase_compare, fermi-lite FFI, BAM error model, and negative tests
make compare-whatshap-phase # optional WhatsHap truth path compared with native phase_compare
```

For private or larger local datasets, pass paths from your shell; no private
paths are embedded in the repository or rendered README.

## README generation

`README.md` is generated from `README.Rmd` so the CLI help and fixture example
outputs stay synchronized with the installed tools:

```bash
make readme
```

This requires R with the `knitr` package. The `readme` target builds the Rust
CLIs first, then renders `README.Rmd`. For a local larger VCF example, set
`PHASE_MNV_EXAMPLE_VCF`, `PHASE_MNV_EXAMPLE_REF`, and
`PHASE_MNV_EXAMPLE_SAMPLE`, then run:

```bash
make readme-external-example
```

## Reference implementations

Reference implementations for comparison/compatibility work are not vendored in
this repository. Clone them when needed. The fermi-lite library is vendored
separately under `vendor/fermi-lite` as an optional local assembly substrate,
not as a phase/MNV reference implementation:

```bash
./scripts/clone_reference_impls.sh
```

### Optional WhatsHap phase-concordance comparison

The external phasing comparison uses WhatsHap plus the native `phase_compare`
binary. It builds a truth path by unphasing the tiny tracked VCF and running
external `whatshap phase`. It builds the query path by running Rust
`phase_mnv_rs --emit all-sites --phase-from-bam` directly on the same VCF/BAM.
Then `phase_compare` reports common heterozygous sites, intersection PS blocks,
assessed adjacent pairs, phase matches, switch errors, and blockwise Hamming
rate:

```bash
make compare-whatshap-phase
```

This is intentionally faster and narrower than hap.py: it measures phasing and
PS/block concordance, not generic variant-call equivalence. The comparison
script strips bcftools/producer command headers and path-bearing
`phase_mnv`/reference records from generated VCFs before comparison, so local
filesystem paths are not retained in comparison VCF headers.

The script discovers `whatshap` from `PATH`, then falls back to a micromamba
environment named `phase-mnv-whatshap`. Override with `WHATSHAP_BIN` or
`WHATSHAP_ENV` when needed.

Detailed semantics and validation notes are in:

- [`docs/semantics.md`](docs/semantics.md)
- [`docs/validation.md`](docs/validation.md)

### Experimental fermi-lite local assembly binding

The repository vendors fermi-lite source under `vendor/fermi-lite` and builds it
through `build.rs` for future local read-evidence adjudication. Bindings are
generated from `fml.h` with bindgen when libclang is available, with a checked-in
narrow fallback for portability (`PHASE_MNV_REQUIRE_BINDGEN=1` forces bindgen).
The fermi-lite wrapper is currently enabled on x86_64 targets because the
vendored upstream `ksw.c` path requires SSE2. The current exposed Rust path is
intentionally small: `src/fermi_lite.rs` wraps
`fml_opt_init`, `fml_assemble`, and
`fml_utg_destroy`, and
`fermi_lite_assemble` is a smoke/utility binary for assembling supplied local
read sequences into FASTA unitigs. It can pass FASTQ/base qualities through to
fermi-lite when `--fastq --ec-k 0` or another non-negative `--ec-k` is used.

`bam_error_model` is a separate helper that learns simple empirical mismatch,
insertion, and deletion summaries from a BAM/CRAM versus a reference, without a
MAPQ filter by default. It reports base-quality and MAPQ bins so a future
`phase_adjudicate` path can separate local assembly from quality-aware evidence
scoring. Its optional `--position-tsv` output is inspired by Heng Li's htsbox
`mapchk`: rows are keyed by one-based read position from the read 5' end
(reverse-strand alignments are reversed) and split into all/low/high/unknown
base-quality groups using `--high-quality-threshold`. Known variant sites are
not masked yet, so mismatch rates from this helper should be interpreted as
sequencing-error-plus-variation unless regions are restricted to
homozygous-reference/high-confidence sites.

Neither helper is yet integrated into `phase_compare` or a final
`phase_adjudicate` workflow.

If you use fermi-lite-backed local assembly results, cite the FermiKit paper
recommended by fermi-lite (Li 2015, Bioinformatics; doi:10.1093/bioinformatics/btv440).

The most relevant upstream tools are:

- `whatshap phase` / `whatshap haplotag` — established read-backed phasing
  reference for BAM/CRAM validation and external phasing workflows.
- `phase_compare` — in-repository native phase-concordance comparator for
  WhatsHap-derived all-sites truth versus Rust BAM-backed all-sites output.
- `vt normalize` — normalization reference for left-aligned + parsimonious
  representation.
- `bcftools norm` — useful validator for emitted normalized records.
- `Illumina Nirvana` — codon/transcript-aware SNV-only MNV recomposition
  reference for future benchmarking of annotation-driven recomposition rules.
- `fermi-lite` / FermiKit — optional local reassembly substrate for future
  read-evidence adjudication of difficult phasing disagreements.

## CI

GitHub Actions builds and tests the Rust project on Linux and macOS:

- Rust release/static where supported
- behavior and negative/failure-mode tests for the Rust binaries
- Rust `bcftools norm` validation of emitted normalized records
- Linux WhatsHap-derived all-sites phase comparison with native `phase_compare`
- native `unphase_vcf`, fermi-lite FFI, and BAM error-model smoke coverage
- binary artifact upload with SHA256 sums

## Notes

- Reads VCF/BCF via htslib/rust-htslib.
- Uses phased diploid GT and `FORMAT/PS`.
- Emits `TYPE=MNV` for pure SNV blocks and `TYPE=COMPLEX` for blocks including indels.
- Normalizes internally with the Tan, Abecasis & Kang 2015 left-aligned + parsimonious rules (doi:10.1093/bioinformatics/btv112).
- Does not require a separate `vt normalize` or `bcftools norm` pass for emitted records.
