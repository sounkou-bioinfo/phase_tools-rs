# phase_tools-rs

`phase_tools-rs` is being reduced to one job:

> **Call medically relevant difficult loci from short-read WGS/WES, emit an
> explicit no-call when the assay cannot support the claim, and attach a
> machine-checkable decision certificate.**

The old public collection of MNV, phasing, contamination, ancestry, BAM error
model, and assembly commands is removed from the build. Useful phasing and
assembly code may return later only as private evidence kernels for a target
caller.

## Closed target scope

The registry is finite. It contains HLA and KIR through the
[Unum](https://github.com/fg-labs/unum) Rust port of T1K, plus every target in
the Illumina DRAGEN v4.5 Targeted Caller set.

| Target | WGS | WES | Backend in this refocus | State |
|---|---:|---:|---|---|
| HLA | yes | capture-dependent | Unum/T1K lane | runnable |
| KIR | yes | capture-dependent | Unum/T1K lane | runnable |
| CYP2B6 | yes | no | native paralog lane | contract only |
| CYP2D6 | yes | no | native paralog lane | contract only |
| CYP21A2 | yes | no | native paralog lane | contract only |
| GBA | yes | no | native paralog lane | contract only |
| HBA | yes | validated enrichment only | native integer solver | solver kernel |
| LPA | yes | no | native repeat-CN lane | contract only |
| RH | yes | no | native blood-group/paralog lane | contract only |
| SMN | yes | validated enrichment only | native paralog lane | contract only |

“Closed” means that target-dependent code must exhaust this enum and the Lean
proof checks the same finite set. It does **not** mean that unfinished callers
are represented as finished. A `contract-only` target cannot issue a valid
`called` certificate.

The WES entries are observability contracts, not marketing labels. HBA and SMN
require a declared, validated target-enrichment profile. Every target still has
read/depth/evidence quality gates after this assay-level check.

## Two algorithmic lanes

### 1. HLA/KIR allele typing

`phase-tools unum` executes the maintained Unum backend without a shell. Unum
owns candidate-read extraction, allele k-mers, banded alignment, abundance
estimation, and allele inference. This repository owns:

- the HLA/KIR-only backend boundary;
- WGS/WES observability declarations;
- input, resource, and result hashes;
- normalized call cardinality;
- the proof-carrying certificate.

This avoids copying the T1K port into a second codebase. Direct `unum-core`
embedding should wait for a stable end-to-end library API; the current
high-level driver lives in the Unum binary crate.

### 2. Copy-number/paralog/repeat targets

The native lane consumes typed evidence rather than pretending all loci share
one pileup caller. Its evidence vocabulary includes unique and total depth,
paralog-differentiating sites, junction reads, small variants, repeat-spanning
reads, and phase links.

The first executable kernel is HBA hypothesis selection. It takes a versioned,
resource-defined hypothesis catalogue and integer-valued evidence. Candidate
penalties are computed as:

```text
prior_penalty
  + sum(ceil(abs(observed - expected) / tolerance) * weight)
```

The unique minimum must beat the runner-up by the requested margin. Otherwise
the result is an explicit no-call. Integer arithmetic makes the exact decision
portable and suitable for certificate verification. Feature extraction,
normalization, population hypothesis resources, and analytical validation are
still separate work; the synthetic example is not a clinical HBA resource.

## Build and test

```bash
make test
make proof
make release
```

Lean is pinned in `lean-toolchain`; Rust has no runtime crate dependencies in
this refocus slice.

## Inspect the target contract

```bash
cargo run -- targets
cargo run -- targets --assay wes
cargo run -- targets --assay wes --validated-enrichment
```

## Run the synthetic HBA solver example

```bash
cargo run -- hba \
  --assay wgs \
  --evidence examples/hba/evidence.synthetic.tsv \
  --hypotheses examples/hba/hypotheses.synthetic.tsv \
  --min-margin 10 \
  --certificate /tmp/hba.cert

cargo run -- verify --certificate /tmp/hba.cert
```

## Run HLA/KIR through Unum

```bash
cargo run -- unum \
  --target HLA \
  --assay wgs \
  --unum /path/to/unum \
  --bam sample.bam \
  --ref-seq hla.ref.fa \
  --ref-coord hla.coord.fa \
  --bam-mode alignment \
  --output-prefix results/sample.hla \
  --threads 8 \
  --certificate results/sample.hla.cert
```

KIR uses the same command with a KIR reference. Resource construction and
versioning remain Unum responsibilities.

## What a Lean certificate proves

The Lean project proves the finite target closure, the DRAGEN-v4.5 targeted
subset, the HBA/SMN WES-enrichment rule, call/no-call cardinality, and soundness
of the HBA winner/margin witness. Rust verifies the corresponding certificate
fields before writing them.

A certificate proves that a result follows the declared deterministic contract
for content-addressed inputs and resources. It cannot prove that the reads came
from the stated patient, that the assay was unbiased, that a resource catalogue
is biologically complete, or that the caller is clinically accurate. Those
claims require sample provenance, truth data, calibration, and validation.

See [the architecture](docs/architecture.md),
[the closed scope](docs/scope.md), and
[the formal assurance boundary](docs/certificates.md).
