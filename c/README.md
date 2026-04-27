# phase_mnv

Minimal phased haplotype builder from a **phased** VCF/BCF plus a FASTA reference.

This deliberately strips away the extra get_MNV functionality (GFF/CDS amino-acid
annotation, BAM read counting, GUI, summaries, etc.). get_MNV is credited as the
motivation for reducing the problem to an MNV-construction primitive; this C
implementation does not copy its Rust code and has a narrower phased-VCF/BCF
scope. It only does the core haplotype operation:

1. read one sample from VCF/BCF with htslib;
2. keep phased diploid alternate alleles (`0|1`, `1|0`, `1|1`, etc.);
3. group variants carried on the same phased haplotype and same `FORMAT/PS`
   phase set;
4. fetch the merged REF span from FASTA;
5. apply the phased alleles to build a constructed VCF record.

Pure SNV blocks are emitted as `TYPE=MNV`; blocks that include an insertion or
deletion are emitted as `TYPE=COMPLEX`. The constructed record is normalized
internally using the Tan, Abecasis & Kang 2015 left-aligned + parsimonious
algorithm (doi:10.1093/bioinformatics/btv112), so no extra `vt normalize` or
`bcftools norm` pass is required for the emitted record.
By default only adjacent records are merged. Use `--max-gap N` if you want to
allow unchanged reference bases between phased variants.

## Build

```bash
make
make install   # installs to ~/.local/bin by default
```

The Makefile uses `pkg-config htslib` if available. On this machine it falls
back to the Rhtslib headers and `/usr/lib/x86_64-linux-gnu/libhts.so.3`.
Override paths if needed:

```bash
make HTS_CPPFLAGS='-I/path/to/htslib/include' HTS_LIBS='-L/path/to/lib -lhts'
```

## Run

```bash
./phase_mnv -r reference.fa -s SAMPLE phased.vcf > mnv.vcf
./phase_mnv -r reference.fa -s SAMPLE phased.bcf -o mnv.vcf
```

Useful options:

```text
-g, --max-gap N      allow up to N unchanged reference bases between variants
    --min-vars N     minimum source variants per emitted call (default 2)
    --min-snvs N     alias for --min-vars
    --unsupported-alleles skip|fail
                      skip selected unsupported alleles, or fail fast
    --warn-on-n      warn when selected REF/ALT contains N
    --no-ref-check   ignore VCF REF vs FASTA mismatches
    --no-header      emit records only
-q, --quiet          suppress stderr summary
```

## Output

The output VCF has one selected sample and these fields:

- `REF` is fetched from the reference across the merged span, then normalized.
- `ALT` is the same span with phased alleles applied, then normalized.
- `INFO/TYPE`: `MNV` for pure SNV blocks, `COMPLEX` when indels are included.
- `INFO/NVAR`: number of original source variants merged.
- `INFO/NSNPS`: number of source SNVs merged.
- `INFO/END`: merged interval end coordinate.
- `INFO/SOURCE_POS`: original source variant positions.
- `INFO/HAPS`: one-based haplotype(s) carrying the merged call (`1`, `2`, or `1,2`).
- `INFO/PS` and `FORMAT/PS`: phase set when present in the input.
- `FORMAT/GT`: constructed phased genotype (`1|0`, `0|1`, `1|1`).

## Scope

Implemented intentionally narrow:

- one selected sample;
- diploid phased GT only;
- plain DNA alleles only (`REF` and selected `ALT` contain A/C/G/T/N; symbolic, breakend, spanning-deletion, and non-DNA selected alleles are skipped by default or rejected with `--unsupported-alleles fail`);
- emits plain text VCF (pipe to `bgzip` if compressed output is needed).

This is the reusable MNV-construction slice; no codon or amino-acid annotation is
attempted here.

## Citation

The normalization implemented here is based on:

> Tan A, Abecasis GR, Kang HM. **Unified representation of genetic variants.**
> *Bioinformatics.* 2015;31(13):2202-2204.
> doi:[10.1093/bioinformatics/btv112](https://doi.org/10.1093/bioinformatics/btv112).

Please cite this paper when relying on the left-aligned + parsimonious
normalized `REF`/`ALT` representation emitted by this tool.
