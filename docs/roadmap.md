# Roadmap

## Gate 0 — scope reset

Complete in this refocus:

- one public binary;
- closed ten-target registry;
- old toolbox binaries removed from the build;
- HLA/KIR Unum adapter;
- HBA integer selection kernel;
- canonical content-addressed certificates;
- Lean registry and certificate proofs.

## Gate 1 — HLA/KIR end-to-end fixtures

- pin Unum and reference-resource versions;
- add public WGS and WES fixtures;
- normalize per-locus genotype and copy-number output;
- define locus-specific evidence/no-call gates;
- differential-test wrapper output while treating upstream output as evidence,
  not truth.

## Gate 2 — HBA extraction

- implement unique and total-depth bins;
- implement panel/cohort normalization with explicit batch metadata;
- collect HBA1/HBA2 differentiating-site likelihoods;
- collect known deletion, duplication, and hybrid junctions;
- emit the prepared integer evidence consumed by the current solver.

## Gate 3 — HBA hypothesis model

- compile chromosome haplotypes into diploid hypotheses;
- add small variants and phase links;
- estimate calibrated integer penalties from truth data;
- version and hash all resources;
- prove the scoring fold and selection implementation against the Lean model.

## Gate 4 — close the native DRAGEN-targeted set

Promote one target at a time from `contract-only`. SMN comes first, followed by
GBA/CYP21A2, CYP2D6/CYP2B6, RH, and LPA. No target is promoted solely because a
generic interface exists.

## Validation gate for every target

- synthetic invariants and negative fixtures;
- public or independently shareable truth where available;
- reproducibility across thread count and input order;
- WGS/WES assay-stratified no-call rates;
- ancestry-aware error analysis;
- explicit comparison to current DRAGEN output without using DRAGEN as the
  sole ground truth.
