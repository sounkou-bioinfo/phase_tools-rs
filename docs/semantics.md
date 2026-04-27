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

## Output model

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

Skipped symbolic/non-DNA alleles are **not currently barriers**. They are ignored
as observations, and nearby supported variants may still merge across their
positions when `--max-gap` permits the unchanged reference span. With
`--max-gap 0`, directly adjacent supported variants are still required.

Tracked regression fixture:

```text
tests/fixtures/symbolic.vcf
tests/fixtures/symbolic.max1.expected.body.vcf
```

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

1. Should unsupported symbolic/breakend/spanning-deletion alleles act as hard
   barriers that prevent merging across their span or position?
2. Should `N` continue to be treated as plain DNA for merging, or should records
   containing ambiguous bases be skipped by default?
3. Should future output support true multi-ALT records when two haplotypes carry
   different ALT haplotypes over the same normalized span?
4. Should overlapping phased records fail loudly instead of silently splitting
   into smaller blocks that may be below `--min-vars`?

Until these are intentionally changed, the tracked fixtures above define the
current behavior.
