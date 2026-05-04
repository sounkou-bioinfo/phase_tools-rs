# phase_tools-rs

`phase_tools-rs` is a Rust genomics toolkit for phased VCF/BCF handling,
read-backed phase evidence, BAM/CRAM empirical error summaries, and
anchor-based contamination probes. The flagship binary, `phase_mnv_rs`, builds
normalized phased haplotype MNV/COMPLEX records from variants on the same
haplotype and phase set.

The code uses `rust-htslib` for VCF/BCF/BAM/CRAM I/O. Rust-only phasing,
adjudication, error-model, contamination, and fermi-lite assembly helpers are
experimental unless a document explicitly says otherwise.

## Tools

| Binary | Purpose |
| --- | --- |
| `phase_mnv_rs` | Build normalized phased `TYPE=MNV` / `TYPE=COMPLEX` records. |
| `phase_compare` | Native phase-concordance and switch-error summary for two VCF/BCF files. |
| `unphase_vcf` | Convert phased VCF/BCF GT separators to unphased and optionally drop phase tags. |
| `bam_error_model` | Summarize empirical mismatch/insertion/deletion patterns from BAM/CRAM vs FASTA. |
| `phase_adjudicate` | Initial read-evidence adjudication for `phase_compare --pair-tsv` SNV pairs. |
| `bam_contamination` | Experimental anchor-site contamination probe from BAM/CRAM plus TSV/VCF/BCF anchors. |
| `bam_ancestry` | Experimental Summix-style ancestry mixture probe from BAM/CRAM plus population-AF anchors. |
| `fermi_lite_assemble` | Small fermi-lite FFI smoke/utility binary for local assembly experiments. |

Full generated CLI help lives in [`docs/cli.md`](docs/cli.md).

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

Useful build targets:

```bash
make test                  # tracked-fixture test suite
make install               # install CLIs to ~/.local/bin
make static-release        # static Linux release when supported
make compare-whatshap-phase # optional Linux WhatsHap comparison
```

Output format is inferred from `-o/--output`: `.vcf`, `.vcf.gz`, `.vcf.bgz`, or
`.bcf`; stdout is plain VCF. Use `--threads N` / `-@ N` for htslib/BGZF worker
threads.

## Minimal MNV example

```bash
target/release/phase_mnv_rs -r tests/fixtures/ref.fa --no-header -q tests/fixtures/phased_mnv.vcf
```

```text
chr1	1	.	AC	GT	.	PASS	TYPE=MNV;NVAR=2;NSNPS=2;END=2;SOURCE_POS=1,2;HAPS=2;PS=10	GT:PS	0|1:10
chr1	4	.	TA	CG	.	PASS	TYPE=MNV;NVAR=2;NSNPS=2;END=5;SOURCE_POS=4,5;HAPS=1;PS=10	GT:PS	1|0:10
chr1	6	.	CG	AT	.	PASS	TYPE=MNV;NVAR=2;NSNPS=2;END=7;SOURCE_POS=6,7;HAPS=1,2;PS=20	GT:PS	1|1:20
```

## BAM-backed examples and helpers

`--phase-from-bam` performs experimental Rust read-backed phasing before MNV
construction. The default `--phase-algorithm mec` solves an exact diploid
single-sample MEC objective per connected component after deterministic read
selection; use WhatsHap itself when established WhatsHap behavior is required.

`bam_error_model` estimates error-plus-variation unless you restrict to trusted
homozygous-reference sites or use its optional site guard:

```bash
target/release/bam_error_model \
  --reference tests/fixtures/ref.fa \
  --position-tsv read_pos.tsv \
  --event-tsv events.tsv \
  tests/fixtures/read_phase.bam
```

`bam_contamination` accepts headered TSV anchors (`chrom,pos,ref,alt,gt[,ref_af]`)
or VCF/VCF.GZ/VCF.BGZ/BCF anchors with sample `GT`. VCF `INFO/REF_AF` is used
as reference allele frequency when present; otherwise `INFO/AF` is treated as
ALT frequency and converted to `ref_af = 1 - AF`. Anchor REF bases are validated
against the supplied FASTA. No MAPQ/baseQ filter is applied by default; optional
thresholds must be explicit.

`bam_ancestry` is an experimental Summix-style helper. It counts REF/ALT bases at
ancestry-informative anchors, estimates observed ALT fractions, and fits a
non-negative mixture over reference population ALT frequencies:

```bash
target/release/bam_ancestry \
  --reference ref.fa \
  --bam reads.bam \
  --anchors ancestry_anchors.tsv \
  --populations AFR,EUR,EAS,SAS,AMR
```

For targeted HLA/HLA-DRB1 exon 2 assays, treat contamination estimates as anchor
probes rather than definitive sample-wide estimates unless the assay includes
independent contamination markers outside the highly polymorphic HLA interval.
See [`docs/contamination_and_ancestry.md`](docs/contamination_and_ancestry.md).

## Citations

If you use this software, cite `phase_tools-rs` via [`CITATION.cff`](CITATION.cff)
and cite the upstream methods relevant to your workflow:

- Normalization / parsimonious REF/ALT representation: Tan, Abecasis & Kang
  2015, *Bioinformatics*, doi:[10.1093/bioinformatics/btv112](https://doi.org/10.1093/bioinformatics/btv112).
- External WhatsHap phasing workflows: Martin et al. 2016, *bioRxiv*,
  doi:[10.1101/085050](https://doi.org/10.1101/085050).
- MEC / weighted haplotype assembly ideas behind the experimental Rust BAM
  phaser: Patterson et al. 2015, *Journal of Computational Biology*,
  doi:[10.1089/cmb.2014.0157](https://doi.org/10.1089/cmb.2014.0157).
- fermi-lite / FermiKit-backed local assembly experiments: Li 2015,
  *Bioinformatics*, doi:[10.1093/bioinformatics/btv440](https://doi.org/10.1093/bioinformatics/btv440).
- CHARR-like contamination concepts: Lu et al. 2023, *American Journal of Human
  Genetics*, doi:[10.1016/j.ajhg.2023.10.011](https://doi.org/10.1016/j.ajhg.2023.10.011).
- Summix-style summary allele-frequency ancestry deconvolution concepts:
  Arriaga-MacKenzie et al. 2021, *American Journal of Human Genetics*,
  doi:[10.1016/j.ajhg.2021.05.016](https://doi.org/10.1016/j.ajhg.2021.05.016).

The Rust BAM phaser and contamination helper are inspired by these methods but
are not drop-in replacements for WhatsHap, CHARR, VerifyBamID, or Summix.

## Development

`README.md` and `docs/cli.md` are generated from R Markdown sources:

```bash
make readme
```

Default tests and CI use only tracked fixtures under `tests/fixtures/`; pass
private/larger paths from your shell rather than committing them.

## License

MIT. The fermi-lite source under `vendor/fermi-lite` carries its upstream
license and is used as an optional local assembly substrate.

## Documentation

- [`docs/cli.md`](docs/cli.md) — generated full CLI help.
- [`docs/semantics.md`](docs/semantics.md) — variant, phasing, and helper semantics.
- [`docs/validation.md`](docs/validation.md) — tracked validation strategy.
- [`docs/contamination_and_ancestry.md`](docs/contamination_and_ancestry.md) — contamination/HLA/ancestry notes.
- [`docs/nirvana_benchmark.md`](docs/nirvana_benchmark.md) — Nirvana-style codon benchmark notes.
