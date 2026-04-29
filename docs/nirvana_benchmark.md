# Nirvana-style MNV recomposition benchmark

Reference implementation under local ignored clone:

```text
reference_impls/Nirvana
HEAD 62d30326985a6de6870d6dd174abeb56474cf466
```

`phase_mnv_rs` currently implements the first narrow benchmark slice as:

```bash
--mnv-algorithm nirvana-codon --codon-map codons.tsv
```

The codon map is BED-like:

```text
CHROM START0 END0 TRANSCRIPT CODON_ID [ignored...]
```

For local GRCh37 experiments, generate an ignored Ensembl GRCh37 codon map from
GTF CDS features:

```bash
./scripts/download_grch37_codon_map.sh
```

For a real VCF, prefer a variant-position subset so startup and memory stay
reasonable:

```bash
VCF=input.vcf.gz ./scripts/download_grch37_codon_map.sh
```

This writes to `resources/` by default and is not committed. The full genome map
can be large. The converter emits multiple intervals with the same
`TRANSCRIPT/CODON_ID` when a codon crosses a splice junction.

Implemented scope:

- SNV-only recomposition.
- Seed rule: at least two phased SNV observations on the same haplotype and
  phase set share a transcript/codon key.
- One pass over the VCF/BCF still collects observations, phase/BAM data, and
  codon memberships; recomposition is then generated from the in-memory fused
  observations.
- Output remains normalized `TYPE=MNV` VCF/BCF records.

Not yet full Nirvana parity:

- no Nirvana transcript-cache reader yet;
- no adjacent-codon aggregation yet;
- no multi-sample/sample-specific recomposition yet;
- homozygous variants are not yet made available to all phase sets;
- unsupported overlapping indel/SV barriers are not yet Nirvana-compatible;
- no `linkedVids`, `isDecomposedVariant`, `isRecomposedVariant`, QUAL/FILTER/GQ
  recomposition compatibility yet.

Tracked fixture:

```text
tests/fixtures/nirvana_codon.vcf
tests/fixtures/nirvana_codon.codons.tsv
tests/fixtures/nirvana_codon.expected.body.vcf
```
