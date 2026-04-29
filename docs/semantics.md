# Variant semantics

This document defines the currently supported semantics for `phase_mnv_rs`. The
C implementation is kept byte-identical for this scope.

## Input model

- Input is VCF/BCF plus FASTA.
- One sample is processed: `--sample NAME`, or the first sample by default.
- Only diploid phased `FORMAT/GT` is used.
  - Accepted examples: `0|1`, `1|0`, `1|1`, `1|2`, `2|1`.
  - Unphased genotypes such as `0/1` are skipped.
  - Missing or non-diploid genotypes are skipped.
- `FORMAT/PS` is honored when present.
  - Variants are only merged within the same haplotype and same phase set.
  - If `PS` is absent/missing, all missing-PS variants share the missing-PS key
    and may merge when phased and close enough.
- By default `--max-gap 0` merges only directly adjacent variants on the same
  haplotype/phase set. `--max-gap N` permits up to `N` unchanged reference bases
  between source variants.

## BAM-backed phasing

The default mode consumes already phased VCF/BCF. The Rust binary also has an
experimental read-backed phasing mode:

```bash
--phase-from-bam reads.bam
```

In this mode input `GT` phase separators and `FORMAT/PS` values are ignored.
The Rust implementation uses `rust-htslib` to read an indexed BAM/CRAM and
extracts per-read allele observations at heterozygous diploid VCF records.

The default phasing algorithm is now:

```bash
--phase-algorithm mec --phase-max-coverage 15
```

This performs deterministic read selection to cap per-variant read coverage,
builds read-connected components, and solves the exact weighted diploid
single-sample MEC objective per component by dynamic programming over active-read
bipartitions. The earlier pairwise parity heuristic remains available as:

```bash
--phase-algorithm greedy
```

This mode is inspired by WhatsHap's non-pedigree read-backed phasing model, but
it is still not a full WhatsHap clone: pedigree/PedMEC, all read-selection
heuristics, all CLI options, and switch-error validated equivalence remain future
work. The compatibility contract is currently the tracked Rust fixture and local
validation against WhatsHap on real data, not byte identity with WhatsHap.

For comparison against the established upstream phaser,
`scripts/phase_from_bam_then_mnv.sh` provides a local workflow that:

1. runs `scripts/unphase_vcf.py` to replace `|` with `/` in GT fields;
2. drops `FORMAT/PS` and `FORMAT/PQ` by default, because those tags describe the
   discarded phase state;
3. runs `whatshap phase` with an explicitly supplied BAM/CRAM and reference
   FASTA;
4. runs `phase_mnv_rs` on the WhatsHap-phased VCF.

The unphasing helper does not change alleles, filters, INFO fields, or non-phase
FORMAT values. It is a preparation step for external phasing comparisons.

## Output model

The default `--emit mnv --mnv-algorithm proximity` mode writes only derived
biallelic merged haplotype records, not the whole input VCF/BCF.

- Output records are biallelic merged haplotype records.
- Pure source-SNV blocks emit `TYPE=MNV`.
- Blocks containing any selected indel or non-SNV plain-DNA allele emit
  `TYPE=COMPLEX`.
- Output includes:
  - `NVAR`: number of source observations in the merged call
  - `NSNPS`: number of source observations that were SNV substitutions
  - `END`: end coordinate of the normalized output REF span
  - `SOURCE_POS`: source VCF positions used in the merged haplotype
  - `HAPS`: haplotypes carrying the same normalized output allele
  - optional `PS`
  - `GT:PS` for the selected sample
- If both haplotypes produce the same normalized merged call, the records are
  merged into one output line with `GT=1|1` and `HAPS=1,2`.
- If haplotypes produce different ALT haplotypes over the same region, they are
  emitted as separate biallelic records, one per ALT haplotype.

The Rust-only `--emit all-sites` mode instead preserves every input VCF/BCF
record and keeps the original input header by duplicating it through htslib, then
appending `phase_mnv` metadata/header records. With `--phase-from-bam` it updates
`FORMAT/GT` and `FORMAT/PS` for phased one-sample inputs; it does not construct
MNV/COMPLEX records in that mode.

## Nirvana-style codon-aware SNV recomposition

The Rust binary has an initial Nirvana-inspired MNV recomposition mode:

```bash
--mnv-algorithm nirvana-codon --codon-map codons.tsv
```

The codon map is a small BED-like text file:

```text
CHROM START0 END0 TRANSCRIPT CODON_ID [ignored...]
```

`START0`/`END0` are 0-based half-open codon spans. A phased SNV observation may
belong to multiple transcript/codon keys. In this first slice, the tool emits an
MNV only when at least `--min-vars` phased SNV observations on the same haplotype
and phase set share a codon key. Indels and non-SNV observations are ignored by
this mode, so emitted calls are `TYPE=MNV` only.

This captures Nirvana's core same-codon SNV seed rule in a lightweight,
validation-friendly form without vendoring Nirvana transcript caches. Full
Nirvana parity still requires transcript-cache integration, adjacent-codon
aggregation, sample-specific multi-sample recomposition, homozygous-variant
availability across phase sets, unsupported-overlap barriers, and exact
`linkedVids`/QUAL/FILTER/GQ behavior.

## Multi-allelic input sites

Multi-allelic input is supported only through the allele selected by the chosen
sample's phased genotype.

For a record like:

```text
REF=A ALT=G,T GT=1|2
```

- haplotype 1 observes `A>G`
- haplotype 2 observes `A>T`
- unselected ALT alleles are ignored
- output remains biallelic; the tool does not emit one multi-ALT output record

This means multi-allelic input is effectively decomposed into haplotype-specific
biallelic observations for the selected sample. Each observation can merge with
nearby observations on the same haplotype and same phase set.

Current caveats:

- A single multi-allelic site alone is not emitted because the default and
  minimum supported `--min-vars` is 2.
- If different haplotypes carry different selected ALT alleles and merge across
  the same positions, the outputs are separate biallelic records unless the
  normalized REF/ALT/position representation is identical.
- Input INFO/FORMAT fields other than `GT`/`PS` are not propagated into the
  output record.

Tracked regression fixture:

```text
tests/fixtures/multiallelic.vcf
tests/fixtures/multiallelic.expected.body.vcf
```

## Symbolic, breakend, spanning-deletion, and non-DNA alleles

The tool only turns plain DNA alleles into haplotype observations. Plain DNA is
currently defined as strings containing only `A`, `C`, `G`, `T`, or `N` after
uppercasing.

Skipped ALT examples:

- symbolic alleles: `<DEL>`, `<INS>`, `<DUP>`, `<NON_REF>`
- spanning deletion allele: `*`
- breakends: `A[chr2:10[`, `]chr2:10]A`
- any allele containing non-DNA characters

If the REF allele itself is not plain DNA, the whole record is skipped.

Default policy is:

```bash
--unsupported-alleles skip
```

Under the default, skipped symbolic/non-DNA alleles are **not currently
barriers**. They are ignored as observations, and nearby supported variants may
still merge across their positions when `--max-gap` permits the unchanged
reference span. With `--max-gap 0`, directly adjacent supported variants are
still required.

Use fail-fast policy to reject selected unsupported alleles instead of skipping
them:

```bash
--unsupported-alleles fail
```

This fails when the selected sample/haplotype points to an unsupported ALT, when
the selected ALT allele index is invalid, when a selected ALT equals REF, or when
REF is unsupported. Unselected multi-allelic ALT alleles remain ignored.

Tracked regression fixture:

```text
tests/fixtures/symbolic.vcf
tests/fixtures/symbolic.max1.expected.body.vcf
```

## gVCF and `<NON_REF>` records

`phase_mnv_rs` can read VCF/BCF files that contain gVCF-style records, but it
currently treats them with the same allele-selection rules as any other VCF:

- homozygous-reference block records such as `ALT=<NON_REF>; GT=0/0` contribute
  no alternate haplotype observation;
- records with a real selected ALT plus an unselected `<NON_REF>` allele can
  still contribute the real selected ALT observation;
- a selected `<NON_REF>` or `<*>` allele is unsupported and is skipped by
  default, or fails under `--unsupported-alleles fail`;
- gVCF `END` block spans are not treated as merge barriers.

For conformance work, prefer a variants-only VCF when available. gVCF input is
acceptable as transport, but full gVCF block semantics are not a recomposition
contract for this tool yet.

## Ambiguous `N` bases

`N` is currently treated as plain DNA, so selected alleles containing `N` can
participate in merged calls. Use:

```bash
--warn-on-n
```

to emit a warning for each selected haplotype observation whose REF or ALT allele
contains `N`. The run summary always reports `observations_with_n`, even when
warnings are disabled. The typo alias `--warm-on-n` is accepted but not shown in
help.

Tracked regression fixture:

```text
tests/fixtures/n_base.vcf
tests/fixtures/n_base.expected.body.vcf
```

## Statistics output

Unless `--quiet` is set, summary statistics are written to stderr and include an
explicit output destination:

```text
phase_mnv: input=... reference=... output=stdout sample=...
phase_mnv: settings max_gap=... min_vars=... unsupported_alleles=... warn_on_n=... no_ref_check=... no_header=... output_format=... threads=... emit=... mnv_algorithm=...
phase_mnv: records=... phased_records=... haplotype_variant_observations=... emitted_calls=...
phase_mnv: skipped no_gt=... non_diploid=... missing_gt=... unphased=... ref_hap_alleles=...
phase_mnv: unsupported ref_non_dna=... alt_out_of_range=... alt_symbolic_or_breakend=... alt_spanning_deletion=... alt_non_dna=... alt_same_as_ref=... unsupported_alt_total=...
phase_mnv: multiallelic_records=... observations_with_n=...
```

When output is written to stdout, the label is exactly `output=stdout`; the
summary itself remains on stderr so it does not corrupt VCF stdout.

## Overlapping records

Overlapping observations on the same haplotype are not merged by the current
windowing rule. In practice they split blocks because a later observation must
start after the current observation's reference end. Normalize/decompose complex
overlapping input before using it when that representation matters.

## Normalization

Merged output records are normalized internally using the left-aligned and
parsimonious representation from:

> Tan A, Abecasis GR, Kang HM. Unified representation of genetic variants.
> Bioinformatics. 2015;31(13):2202-2204. doi:10.1093/bioinformatics/btv112.

No external `vt normalize` or `bcftools norm` post-pass is required for the
supported output scope, though both remain useful validators.

## Open semantic decisions

The main unresolved policy questions are:

1. Should unsupported symbolic/breakend/spanning-deletion alleles gain a third
   `barrier` policy that prevents merging across their span or position?
2. Should `N` continue to be treated as plain DNA for merging, or should records
   containing ambiguous bases be skipped or failed by default?
3. Should future output support true multi-ALT records when two haplotypes carry
   different ALT haplotypes over the same normalized span?
4. Should overlapping phased records fail loudly instead of silently splitting
   into smaller blocks that may be below `--min-vars`?

Until these are intentionally changed, the tracked fixtures above define the
current behavior.
