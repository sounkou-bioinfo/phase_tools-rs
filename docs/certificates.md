# Formal assurance and decision certificates

## Certificate format

Certificates are canonical `key=value` records. Version 1 binds:

- target and assay declaration;
- target-enrichment declaration;
- registered backend and backend version;
- called/no-call state and call cardinality;
- SHA-256 of the input, resource bundle, and normalized output;
- an optional HBA selection witness.

The HBA witness contains the candidate count, winner index, winner penalty,
runner-up penalty, and required margin. The verifier accepts it only when:

```text
candidate_count > 0
winner_index < candidate_count
winner_penalty + required_margin <= runner_up_penalty
```

A called HBA result must have exactly one call and a valid witness.

## Exact decision-state contract

A certificate is rejected unless all of the following agree:

- the target and its registered backend;
- WGS/WES observability and the enrichment declaration;
- the target's implementation state;
- called versus no-call cardinality;
- the no-call reason permitted by the assay or implementation state.

Assay failures are exact: an unsupported WES target must report
`assay-not-observable`, while HBA or SMN WES without validated enrichment must
report `missing-validated-enrichment`. These reasons cannot be used for an
observable assay.

For an observable assay, the runnable Unum HLA/KIR lane may report
`insufficient-evidence`. The HBA solver kernel may report either
`insufficient-evidence` or `ambiguous-top-score`. A contract-only target cannot
emit a called result or an evidence-level no-call. It may only participate in a
registry-level assay no-call, such as unsupported WES.

Backend crashes, unreadable inputs, and malformed or mismatched resources are
errors. They are deliberately not encoded as biological no-call reasons.

## Lean project

The proof root is `PhaseTools.lean`.

`PhaseTools.Registry` proves:

- every target constructor occurs in the closed registry list;
- every DRAGEN v4.5 targeted-caller target is inside that registry;
- the DRAGEN targeted WES subset is exactly HBA or SMN;
- only HLA, KIR, and HBA are currently implemented enough to issue calls;
- Unum is registered only for HLA and KIR;
- WGS observability and the HBA/SMN enrichment transition.

`PhaseTools.Certificate` defines the same target, assay, backend,
implementation, call-cardinality, and no-call-reason invariant as the Rust
verifier. It proves that successful verification implies:

- the backend matches the target;
- a called target is implemented;
- a called target is observable under the declared assay;
- a no-call has zero calls;
- an HBA selection winner is in range and satisfies the requested margin.

CI builds with warnings as errors, runs Lean's environment checker, and audits
the compiled `PhaseTools` namespace for axioms outside the standard allowlist.
The audit catches `sorry`, `admit`, `native_decide`, and home-rolled axioms.

## Trusted computing base

The current trusted boundary includes:

- correspondence between the Rust and Lean definitions;
- Rust parsing and SHA-256 implementation;
- the operating system and file reads;
- the target resource creator;
- the evidence extractor or Unum backend;
- assay/sample provenance.

The proofs do not establish biological completeness or clinical validity.
They prevent a narrower but important class of errors: silently calling outside
the assay contract, claiming an unimplemented target, accepting a wrong backend
or no-call reason, accepting an out-of-range winner, accepting an insufficient
selection margin, or confusing zero calls with a called result.

Future hardening should generate Rust and Lean target definitions from one
small declarative source and formally verify the integer HBA scoring fold, not
only its emitted winner witness.
