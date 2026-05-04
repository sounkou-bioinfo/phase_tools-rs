# Contamination and ancestry notes

This repository is being renamed to `phase_tools-rs` because the toolset now
extends beyond MNV construction into read evidence, empirical error modeling,
contamination probes, and future ancestry-aware calibration.

## HLA typing context

For targeted HLA typing, especially HLA-DRB1 exon 2 / shared-epitope analyses,
contamination and sample mix-up checks are important but hard to estimate from
that locus alone:

- HLA-DRB1 exon 2 is highly polymorphic and has strong homology to other HLA
  genes, so apparent low-level alleles can reflect contamination, mapping
  ambiguity, paralogous sequence, allele dropout, or true allelic imbalance.
- A narrow HLA target usually contains too few independent sites to estimate
  genome-wide contamination or ancestry robustly.
- A useful QC package should therefore combine HLA-specific coverage/alignment
  review with independent anchor sites when available.

## Contamination estimation strategy

Empirical sequencing-error rates and contamination are related but not
identical. Estimating an error model from a contaminated BAM can absorb some
contamination signal into the apparent error rate, even when obviously variable
sites are skipped. For contamination, the stronger approach is to use anchor
sites with known expected genotype states or population frequencies.

Practical tiers:

1. **Variant-call based CHARR-like signal**: at high-confidence homozygous-alt
   variants, count reference read infiltration using AD-like data. This is cheap
   when VCF/gVCF records contain `GT`, `AD`, `DP`, and ideally population
   reference allele frequency.
2. **BAM anchor-site signal**: when BAM/CRAM is available, count REF/ALT bases
   at a curated anchor panel. This is what the experimental `bam_contamination`
   helper starts to do.
3. **Full model-based tools**: for production genome/exome contamination, use
   established tools such as VerifyBamID2/ContEst/Conpair where appropriate,
   especially when population frequencies or ancestry are uncertain.

For targeted HLA-only data, `bam_contamination` should be interpreted as an
anchor-site probe rather than a definitive genome-wide contamination estimator.
It becomes more useful when the assay includes independent SNP anchors outside
HLA or a well-designed target-panel contamination marker set.

## Current `bam_contamination` scope

`bam_contamination` consumes:

```text
--bam reads.bam|reads.cram
--reference ref.fa
--anchors anchors.tsv|anchors.vcf|anchors.vcf.gz|anchors.vcf.bgz|anchors.bcf
```

Anchor TSV input must include a header. Required columns:

```text
chrom  pos  ref  alt  gt  [ref_af]
```

Anchor VCF/VCF.GZ/VCF.BGZ/BCF input is also supported. The selected sample's `FORMAT/GT`
defines the expected genotype; `INFO/REF_AF` is interpreted as reference allele
frequency when present, and `INFO/AF` is interpreted as alternate allele
frequency otherwise (`ref_af = 1 - AF`). Only biallelic SNV anchors are used.
REF alleles from both TSV and VCF/BCF anchors are validated against the supplied
FASTA so coordinate/build mismatches fail fast.

It reports per-anchor REF/ALT/other counts and raw reference infiltration at
homozygous-alt anchors. If `ref_af` is supplied, it also reports a simple
CHARR-like component:

```text
(ref_reads / (ref_reads + alt_reads)) / ref_af
```

Per-anchor `ref_fraction` and `alt_fraction` use `ref_reads + alt_reads` as the
denominator. `other_fraction` uses `ref_reads + alt_reads + other_reads` so it
captures non-allelic read support separately from REF/ALT balance.

This is intentionally labelled CHARR-like, not a complete CHARR implementation.
Full CHARR uses carefully filtered autosomal variants and population-frequency
adjustment over many sites.

## Ancestry inference plan

Summix-style ancestry deconvolution is a good fit for `phase_tools-rs`, but it
is a separate problem from sequencing contamination. Summix estimates mixture
proportions by minimizing the squared difference between observed allele
frequencies and a convex combination of reference-population allele frequencies:

```text
minimize || AF_observed - R * pi ||^2
subject to pi_k >= 0 and sum(pi_k) = 1
```

A Rust implementation can start with projected-gradient or active-set quadratic
optimization over a small number of reference ancestries. Inputs should be
mainstream indexed tables first:

```text
ref_freqs.tsv.bgz + .tbi
```

Suggested columns:

```text
chrom  pos  ref  alt  AFR  EUR  EAS  SAS  AMR/IAM ...
```

Current and possible tools:

- `bam_ancestry`: experimental BAM/CRAM ancestry mixture probe. It counts
  REF/ALT bases at ancestry-informative SNV anchors and fits non-negative
  reference-population mixture proportions from observed ALT fractions. Anchor
  TSV columns are `chrom pos ref alt POP1 POP2 ...`; population columns are ALT
  frequencies and can be selected with `--populations`.
- `ancestry_mix`: future summary-table estimator for sample/panel ancestry
  proportions from observed AFs and tabix-indexed reference frequencies.
- `ancestry_adjust_af`: future helper to apply Summix-style ancestry adjustment
  to external control AFs for a requested target ancestry mixture.

For HLA typing, genome-wide ancestry estimates are often more stable than local
HLA-only estimates. Local ancestry around HLA may differ from genome-wide
ancestry due to selection and admixture, so any HLA-local ancestry inference
should be reported as exploratory unless supported by a dense, validated panel.

## Indexed reference resources

Use interoperable formats first:

- BGZF + tabix/CSI for reference-frequency and anchor tables.
- BED-like BGZF + tabix/CSI for target intervals.
- VCF/BCF when allele/genotype semantics matter.

A future `cgranges`-style in-memory interval layer can accelerate repeated
overlap checks after loading indexed interval chunks, but the source of truth
should remain portable BGZF/tabix tables until profiling shows a stronger cache
is needed.

## Privacy notes

Contamination and ancestry outputs can reveal sensitive genetic information,
especially on small targeted panels. Avoid embedding private paths, sample names,
or raw per-read data in committed examples. Aggregate outputs should still be
handled as potentially identifying when the anchor panel is small or clinically
sensitive.
