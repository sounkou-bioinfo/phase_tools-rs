# Architecture

The repository has one public binary, `phase-tools`, and one closed target
registry. The design rejects the former “one more utility” growth pattern.

```text
BAM/CRAM + assay declaration + versioned target resources
                         |
                 target-specific evidence
                         |
       +-----------------+------------------+
       |                                    |
Unum allele-typing lane           native target lane
(HLA, KIR)                        (HBA, then DRAGEN set)
       |                                    |
       +--------------- normalized call ----+
                         |
        call / typed no-call + content hashes
                         |
              decision certificate verifier
                         |
                    Lean contract
```

## Public modules

- `model`: closed target, assay, backend, evidence, and no-call types.
- `registry`: exhaustive target contracts and observability checks.
- `unum`: shell-free HLA/KIR backend adapter.
- `hba`: resource-driven integer hypothesis solver.
- `digest`: dependency-free SHA-256 and named-resource hashing.
- `certificate`: canonical certificate format and verifier.

Phasing, local assembly, pileup extraction, and read realignment are
implementation details. They become public only if an external contract
requires them; otherwise they remain target-lane kernels.

## Target implementation rule

A target progresses through three explicit states:

1. `contract-only`: its assay/evidence/output contract exists, but it cannot
   issue a called certificate.
2. `solver-kernel`: the deterministic inference kernel exists for prepared
   evidence.
3. `runnable`: extraction, inference, output normalization, certificates, and
   validation fixtures are wired end to end.

Changing a status requires tests and a registry/proof review. There is no
generic “experimental caller” state that silently produces plausible output.

## Evidence before algorithms

Each target defines the evidence it can legitimately consume. Shared
algorithms are extracted only after two target implementations demonstrate the
same semantics. This avoids prematurely forcing HLA allele abundance, HBA copy
number, LPA repeats, and RH hybrids through one generic caller abstraction.

## Output contract

Every attempted target returns either:

- a call with nonzero call cardinality; or
- a typed no-call.

A software crash, malformed resource, or unreadable input is an error, not a
biological no-call. Assay-not-observable and missing-enrichment results are
no-calls because they are properties of the declared experiment.

Every emitted certificate binds the input, resource set, and normalized output
with SHA-256.
