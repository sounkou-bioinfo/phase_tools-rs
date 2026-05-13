# phase_tools-rs CLI help

This file is generated from `docs/cli.Rmd` by `make readme`. It mirrors the
`--help` text from the current release binaries.

## `phase_mnv_rs`

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
      --write-index[=FMT]
                        Build an index after writing -o output. Indexable
                        -o outputs self-index by default with CSI; FMT can
                        be csi or tbi (BGZF VCF only)
      --no-write-index  Do not build an index for indexable -o output
      --emit MODE        Output mode: mnv (default), combined, or all-sites.
                        combined emits merged MNV/COMPLEX records plus
                        selected-sample input variants not represented by a
                        merge; all-sites preserves input records/header and
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
                        BAM phasing algorithm: mec or greedy (default: mec;
                        whatshap is an alias for mec)
      --tag TAG         Phasing tag for --emit all-sites and constructed output:
                        PS (default) or HP
      --only-snvs       Phase only biallelic SNV genotypes from BAM/CRAM
      --output-read-list FILE
                        Write selected BAM/CRAM reads used for MEC phasing
      --ignore-read-groups
                        Ignore BAM read-group sample names and use all reads
      --use-supplementary
                        Include supplementary alignments in BAM phasing
      --supplementary-distance N
                        Accepted for WhatsHap CLI compatibility; read grouping
                        is currently by QNAME regardless of this distance
      --phase-realign-overhang N
                        REF/ALT realignment overhang for BAM allele detection
                        (default: 10, matching WhatsHap reference mode)
      --phase-max-coverage N
                        Maximum selected read coverage per variant span for MEC
                        phasing (default: 15; alias: --phase-internal-downsampling;
                        inspired by WhatsHap --internal-downsampling)
      --phase-min-mapq N  Minimum read MAPQ for --phase-from-bam (default: 20;
                        aliases: --mapping-quality, --mapq)
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
    exact diploid single-sample MEC per connected component after WhatsHap-style
    QNAME read-pair grouping, deterministic read selection, interval coverage
    capping, bridge rescue, and blank-aware active-read DP; greedy keeps the
    earlier pairwise parity heuristic. Phase comparisons should use switch or
    pairwise phase metrics because a whole-block 0|1/1|0 orientation flip is
    equivalent.
  * With the default --max-gap 0 and --mnv-algorithm proximity, only
    adjacent phased variants are merged. Pure SNV blocks are TYPE=MNV;
    blocks containing indels are TYPE=COMPLEX. nirvana-codon mode only
    recomposes phased SNVs sharing a codon key from --codon-map.
  * Output format is inferred from -o/--output. BCF output always includes
    a VCF/BCF header even if --no-header is set.
  * --emit combined keeps the original VCF/BCF header for the selected
    sample, appends phase_mnv metadata, and writes constructed MNV/COMPLEX
    records plus residual selected-sample input records. Partial
    multi-allelic residuals are allele-aware; multi-sample inputs are
    subset to the selected sample, with common cohort INFO tags stripped
    after subsetting. --phase-from-bam is not supported in combined mode.
  * --emit all-sites keeps the original VCF/BCF header via htslib and
    appends phase_mnv metadata instead of replacing it.
  * Indexable -o outputs (.vcf.gz, .vcf.bgz, .bcf) build a CSI sidecar
    by default after the output file is closed. Use --write-index=tbi for
    a tabix/TBI index on BGZF VCF, or --no-write-index to disable.
    Indexing requires coordinate-sorted .vcf.gz/.vcf.bgz/.bcf output.
  * Unless --quiet is set, summary stats go to stderr and include
    input/reference/output (output=stdout for VCF stdout), settings,
    skip counts, unsupported categories, and N counts.
```

## `phase_compare`

```text
usage: phase_compare [options] truth.vcf|bcf query.vcf|bcf

Fast phase-concordance comparison for two VCF/BCF files. The tool compares
exact shared variant records, diploid heterozygous GT/HP phase, PS/HP block
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

## `unphase_vcf`

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

## `fermi_lite_assemble`

```text
usage: fermi_lite_assemble [options] [--seq SEQ ...]

fermi-lite FFI local assembly utility. With --seq, assembles the supplied
sequences. Without --seq, reads one plain sequence per non-empty stdin line,
ignoring FASTA-style header lines. With --fastq, reads FASTQ from stdin and
passes base qualities to fermi-lite's error-correction path when --ec-k >= 0.
This is intended for assembly-backed local adjudication workflows, not as a full
fermi-lite CLI replacement.

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

## `bam_error_model`

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
--event-tsv FILE      Write exact-Q event/ref/read-base composition TSV
--high-quality-threshold N  BaseQ threshold for high/low position groups (default: 20)
--skip-high-nonref-fraction F  Skip likely variant sites above high-Q non-ref fraction F (off)
--include-duplicates  Include duplicate reads
--include-secondary   Include secondary alignments
--include-supplementary Include supplementary alignments
-h, --help                Show this help
```

## `phase_adjudicate`

```text
usage: phase_adjudicate --reference ref.fa --bam reads.bam --variants variants.vcf|bcf --pair-tsv pairs.tsv [options]

Experimental read-evidence adjudicator for phase_compare pair TSV rows. The
first implementation is deliberately narrow: it adjudicates biallelic SNV pairs
by counting reads spanning both sites and comparing observed read allele parity
with the truth/query phased GT patterns from the pair TSV. No MAPQ or baseQ
filter is applied by default; optional thresholds are explicit.

options:
-r, --reference FILE       Required FASTA reference (used for CRAM decoding)
--bam FILE             Indexed BAM/CRAM read evidence
--variants FILE        VCF/BCF containing the pair sites and alleles
--pair-tsv FILE        phase_compare --pair-tsv output
-o, --output FILE          Output TSV (default: stdout)
-@, --threads N            htslib reader threads (default: 1)
--min-mapq N           Optional MAPQ cutoff (default: 0; no cutoff)
--min-baseq N          Optional per-site baseQ cutoff (default: 0; no cutoff)
--include-duplicates   Include duplicate reads
--include-secondary    Include secondary alignments
--include-supplementary Include supplementary alignments
--assembly-fasta FILE  Experimental fermi-lite local assembly sidecar FASTA
--assembly-tsv FILE    Experimental per-unitig assembly parity sidecar TSV
--use-assembly-decision Use assembly evidence to break otherwise ambiguous decisions
--assembly-window N    Bases of padding around each pair for assembly (default: 100)
--assembly-context N   Bases around pair used for unitig parity scoring (default: 10)
--assembly-max-reads N Maximum reads per pair assembly (default: 200)
--assembly-min-asm-ovlp N fermi-lite minimum assembly overlap (default: 21)
-h, --help                 Show this help
```

## `bam_contamination`

```text
usage: bam_contamination --reference ref.fa --bam reads.bam --anchors anchors.tsv|vcf|vcf.gz|vcf.bgz|bcf [options]

Experimental anchor-site contamination probe for BAM/CRAM data. It counts read
bases at caller-supplied anchor sites and reports raw reference infiltration at
homozygous-alternate anchors, with an optional CHARR-like adjustment when a
reference allele frequency is supplied. It applies no MAPQ/baseQ filter by
default; optional thresholds are explicit.

Anchor TSV requires a header with columns: chrom, pos, ref, alt, gt, optional
ref_af. Positions are 1-based. Anchor VCF/VCF.GZ/VCF.BGZ/BCF input is also supported; GT is
read from the selected sample, INFO/REF_AF is used when present, and INFO/AF is
interpreted as ALT frequency when REF_AF is absent. GT currently supports
biallelic 0/0, 0/1, 1/0, and 1/1 forms with / or |. REF alleles are validated
against the supplied FASTA.

options:
-r, --reference FILE       Required FASTA reference (REF validation; CRAM decoding)
--bam FILE             Indexed BAM/CRAM read evidence
--anchors FILE         Anchor TSV or VCF/VCF.GZ/VCF.BGZ/BCF SNV anchors
--sample NAME          Sample name for VCF/BCF anchors (default: first sample)
-o, --output FILE          Output TSV (default: stdout)
-@, --threads N            htslib reader threads (default: 1)
--min-mapq N           Optional MAPQ cutoff (default: 0; no cutoff)
--min-baseq N          Optional baseQ cutoff (default: 0; no cutoff)
--include-duplicates   Include duplicate reads
--include-secondary    Include secondary alignments
--include-supplementary Include supplementary alignments
-h, --help                 Show this help
```

## `bam_ancestry`

```text
usage: bam_ancestry --reference ref.fa --bam reads.bam --anchors ancestry.tsv [options]

Experimental BAM/CRAM ancestry mixture probe. It counts REF/ALT bases at
caller-supplied ancestry-informative SNV anchors, estimates observed ALT
fractions, and fits a constrained least-squares mixture over reference
population ALT frequencies. It applies no MAPQ/baseQ filter by default;
optional thresholds are explicit. This is Summix-style in spirit, but it is not
a full Summix replacement.

Anchor TSV requires a header with columns: chrom, pos, ref, alt, then reference
population ALT-frequency columns. Positions are 1-based. Use --populations to
select/order population columns; otherwise all columns after alt are used. REF
alleles are validated against the supplied FASTA.

options:
-r, --reference FILE       Required FASTA reference (REF validation; CRAM decoding)
--bam FILE             Indexed BAM/CRAM read evidence
--anchors FILE         Anchor TSV with chrom,pos,ref,alt,popAF...
--populations LIST     Comma-separated population columns to use
-o, --output FILE          Output TSV (default: stdout)
-@, --threads N            htslib reader threads (default: 1)
--min-mapq N           Optional MAPQ cutoff (default: 0; no cutoff)
--min-baseq N          Optional baseQ cutoff (default: 0; no cutoff)
--min-observations N   Minimum REF+ALT observations for fitting (default: 1)
--include-duplicates   Include duplicate reads
--include-secondary    Include secondary alignments
--include-supplementary Include supplementary alignments
-h, --help                 Show this help
```
