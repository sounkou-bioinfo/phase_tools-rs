# phase_tools-rs roadmap

This roadmap keeps the implementation staged and validation-driven. The target is
not raw haplotype-label identity: a whole-block flip such as `0|1,0|1` versus
`1|0,1|0` is equivalent. Phasing parity is assessed with switch errors,
pairwise phase, blockwise Hamming rate, and assessed-pair counts.

## Stage 1: WhatsHap-like read-backed phasing

Goal: match WhatsHap's non-pedigree read-backed behavior as closely as possible
for the supported scope, while failing loudly for unsupported modes.

Current supported scope:

- indexed BAM/CRAM phase input
- diploid single-sample heterozygous VCF/BCF records
- MEC-like exact DP per read-connected component
- QNAME read-pair grouping
- internal downsampling aliases (`--phase-max-coverage`,
  `--internal-downsampling`, `--phase-internal-downsampling`)
- interval/span coverage capping
- bridge-read rescue
- blank-aware active-read DP across unobserved intermediate variants
- `PS` and `HP` output tags
- orientation-aware comparison through `phase_compare`

Remaining work before claiming closer WhatsHap parity:

1. Validate the initial padded REF/ALT realignment path on indel and complex
   allele phasing against WhatsHap on synthetic and public truth-style data.
2. Decide whether affine-gap/k-mer allele detection modes are needed, or keep the
   current unit edit-distance path as the supported reference-mode approximation.
3. Add optional read-list details that are closer to WhatsHap's output columns,
   while keeping the current compact TSV stable or versioned.
4. Implement VCF/BCF phase-input superreads, or explicitly keep them out of
   scope.
5. Decide whether to implement WhatsHap ReadMerger or keep only QNAME fragment
   grouping.
6. Keep pedigree/PedMEC, HapChat, and genotype-distrust options as explicit
   unsupported failures until a dedicated implementation exists.

Validation gate:

- tracked fixtures pass
- public or shareable read-backed comparisons show zero or explainable switch
  differences under orientation-aware metrics
- real replicate checks improve or at least do not regress assessed-pair and
  switch-rate summaries

## Stage 2: local assembly as phasing/adjudication evidence

Goal: use local assembly to adjudicate ambiguous read-backed phase and complex
alleles, not to replace the validated phasing path prematurely.

Implementation direction:

1. Reuse the existing fermi-lite FFI as the first local-assembly substrate.
2. Start with region windows and read subsets emitted by the phasing/adjudication
   engine.
3. Emit assembly contigs plus alignment/allele-support diagnostics as sidecars.
   Initial `phase_adjudicate --assembly-fasta` support writes fermi-lite unitigs,
   and `--assembly-tsv` reports simple edit-distance SNV-pair parity scoring for
   manual inspection.
4. Use assembly evidence to resolve difficult MNV/COMPLEX loci only after
   read-backed phase remains ambiguous or conflicting. The first guarded step is
   `phase_adjudicate --use-assembly-decision`, which can break otherwise
   ambiguous SNV-pair decisions only when the assembly TSV sidecar is also
   emitted.
5. Keep fermi-lite x86_64-only constraints explicit.

Validation gate:

- synthetic local assemblies recover expected haplotypes
- assembly adjudication changes are explainable in sidecar evidence
- no default path changes until assembly-assisted decisions improve real-data
  switch/adjudication metrics

## Stage 3: multi-region BAM variant calling

Goal: move from phasing known VCF variants to generating candidate calls from BAM
across many regions, with a reusable library-first core.

Posterior haplotype-caller model direction:

- generate candidate variants per active region
- build candidate haplotypes
- evaluate read likelihoods against haplotypes
- compute genotype/haplotype posterior probabilities
- derive phase sets from posterior-supported relationships
- emit VCF/BCF plus diagnostics

DRAGEN targeted-caller design lessons to incorporate where appropriate:

- specialized treatment for high-homology and repetitive regions
- explicit event grouping (`EVENT`, `EVENTTYPE`) for gene conversions, VNTRs, and
  ambiguous homologous calls
- joint genotypes/copy-number-aware fields for duplicated regions when needed
- distinct filters for low-confidence targeted calls versus ordinary germline
  calls
- clear warnings that targeted/high-homology callers require appropriate WGS-like
  coverage and are not reliable on arbitrary low-coverage or exome-only data

Initial implemented foundation:

- `phase_tools::mrjd` and `multi_region_joint_detect` provide a manifest-driven
  SNV evidence scanner that groups candidate observations by homologous offset
  across user-specified region groups and emits audit TSV diagnostics plus an
  optional diagnostic VCF sidecar with one record per alt-positive region
  observation and `EVENT`/`EVENTTYPE` tags linking homologous evidence. VCF INFO
  string values are percent-escaped. This is a bootstrap candidate-evidence
  layer, not a posterior caller.

Initial non-goals:

- claiming DRAGEN equivalence
- claiming equivalence to any existing posterior haplotype caller
- full targeted pharmacogene/star-allele calling
- high-homology gene conversion calling without target-specific resources and
  validation

Validation gate:

- region-level VCF correctness on synthetic truth
- public truth-set comparisons for ordinary variant calling
- target-specific validation only when target resources, expected outputs, and
  quality thresholds are explicitly defined

## Engineering rules

- Keep private data and paths under ignored `local_runs/`.
- Implement reusable kernels in the `phase_tools` library crate before exposing
  new binary-only workflows.
- Prefer `rust-htslib` and existing vendored fermi-lite over new heavy
  dependencies.
- Keep every claimed compatibility contract tied to commands, versions, and
  validation outputs.
- Add generic manifest-driven scripts only when they do not embed private paths or
  sample names.
