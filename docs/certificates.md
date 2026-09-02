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

## Lean project

The proof root is `PhaseTools.lean`.

`PhaseTools.Registry` proves:

- every target constructor occurs in the closed registry list;
- every DRAGEN v4.5 targeted-caller target is inside that registry;
- the DRAGEN targeted WES subset is exactly HBA or SMN;
- WGS observability and the HBA/SMN enrichment transition.

`PhaseTools.Certificate` proves soundness of the Boolean call-cardinality and
selection-witness verifiers. The Rust verifier mirrors these propositions and
names the Lean contract in every certificate.

CI runs `lake build` with the pinned Lean release.

## Trusted computing base

The current trusted boundary includes:

- correspondence between the Rust and Lean definitions;
- Rust parsing and SHA-256 implementation;
- the operating system and file reads;
- the target resource creator;
- the evidence extractor or Unum backend;
- assay/sample provenance.

The proofs do not establish biological completeness or clinical validity.
They prevent a narrower class of errors: silently calling outside the assay
contract, claiming an unimplemented target, accepting an out-of-range winner,
accepting an insufficient selection margin, or confusing zero calls with a
called result.

Future hardening should generate Rust and Lean target definitions from one
small declarative source and formally verify the integer HBA scoring fold, not
only its emitted winner witness.
