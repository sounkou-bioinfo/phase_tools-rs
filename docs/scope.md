# Closed target scope

## Registry

The code is closed over ten targets:

```text
HLA KIR CYP2B6 CYP2D6 CYP21A2 GBA HBA LPA RH SMN
```

The last eight are the DRAGEN v4.5 Targeted Caller names. HLA is handled by a
separate allele-typing lane; KIR is added beside it because Unum/T1K supports
the same reference/abundance strategy.

This project uses DRAGEN as a target and output-behaviour comparator, not as a
source of truth and not as an architecture to clone blindly.

## Assay observability

All ten targets may be attempted from suitable short-read WGS after ordinary
sample and locus quality checks.

For WES:

- HLA and KIR may be attempted when the capture supplies the required allele
  evidence; evidence-level gates still decide call versus no-call.
- HBA and SMN require a validated target-enrichment declaration.
- CYP2B6, CYP2D6, CYP21A2, GBA, LPA, and RH are outside the current WES
  contract.

The registry intentionally distinguishes assay observability from caller
implementation status. A target can be observable while its native caller is
still `contract-only`.

## Target slices

### HLA and KIR

Backend: Unum.

Required work owned here: reference manifest, invocation, normalized output,
quality policy, call cardinality, hashes, and certificate. Allele inference
stays upstream in Unum.

### HBA

Current slice: deterministic hypothesis selection from prepared integer
features.

Required next work:

1. unique/total-depth feature extraction with cohort or panel normalization;
2. HBA1/HBA2 differentiating-site likelihoods;
3. deletion/duplication/hybrid junction evidence;
4. small-variant evidence and phase links;
5. a population-complete, versioned chromosome-haplotype catalogue;
6. WGS and enriched-WES calibration;
7. independent truth validation.

### Remaining DRAGEN-targeted loci

Implementation order should follow evidence reuse, not alphabetical order:

1. SMN, because it shares total-copy/paralog-site structure with HBA;
2. GBA and CYP21A2, because they exercise homolog/pseudogene conversion;
3. CYP2D6 and CYP2B6, adding star-allele and hybrid structure;
4. RH, adding blood-group haplotypes and hybrid genes;
5. LPA, adding repeat copy number.

Each target must reach a target-specific validation gate before its registry
state changes from `contract-only`.
