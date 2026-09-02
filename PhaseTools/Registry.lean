namespace PhaseTools

inductive Target where
  | hla
  | kir
  | cyp2b6
  | cyp2d6
  | cyp21a2
  | gba
  | hba
  | lpa
  | rh
  | smn
deriving DecidableEq, Repr

inductive Assay where
  | wgs
  | wes
deriving DecidableEq, Repr

inductive WesSupport where
  | standardCapture
  | validatedTargetEnrichment
  | unsupported
deriving DecidableEq, Repr

def allTargets : List Target :=
  [.hla, .kir, .cyp2b6, .cyp2d6, .cyp21a2, .gba, .hba, .lpa, .rh, .smn]

def dragen45Targets : List Target :=
  [.cyp2b6, .cyp2d6, .cyp21a2, .gba, .hba, .lpa, .rh, .smn]

def dragen45WesTargets : List Target :=
  [.hba, .smn]

def wesSupport : Target → WesSupport
  | .hla => .standardCapture
  | .kir => .standardCapture
  | .hba => .validatedTargetEnrichment
  | .smn => .validatedTargetEnrichment
  | .cyp2b6 => .unsupported
  | .cyp2d6 => .unsupported
  | .cyp21a2 => .unsupported
  | .gba => .unsupported
  | .lpa => .unsupported
  | .rh => .unsupported

def observable (target : Target) (assay : Assay) (validatedEnrichment : Bool) : Bool :=
  match assay with
  | .wgs => true
  | .wes =>
      match wesSupport target with
      | .standardCapture => true
      | .validatedTargetEnrichment => validatedEnrichment
      | .unsupported => false

theorem mem_allTargets (target : Target) : target ∈ allTargets := by
  cases target <;> simp [allTargets]

theorem dragen45_closed (target : Target) (_membership : target ∈ dragen45Targets) :
    target ∈ allTargets := by
  exact mem_allTargets target

theorem dragen45_wes_iff (target : Target) :
    target ∈ dragen45WesTargets ↔ target = .hba ∨ target = .smn := by
  cases target <;> simp [dragen45WesTargets]

theorem all_wgs_observable (target : Target) :
    observable target .wgs false = true := by
  rfl

theorem hba_plain_wes_not_observable :
    observable .hba .wes false = false := by
  rfl

theorem hba_enriched_wes_observable :
    observable .hba .wes true = true := by
  rfl

theorem smn_plain_wes_not_observable :
    observable .smn .wes false = false := by
  rfl

theorem smn_enriched_wes_observable :
    observable .smn .wes true = true := by
  rfl

theorem hla_wes_observable :
    observable .hla .wes false = true := by
  rfl

theorem kir_wes_observable :
    observable .kir .wes false = true := by
  rfl

end PhaseTools
