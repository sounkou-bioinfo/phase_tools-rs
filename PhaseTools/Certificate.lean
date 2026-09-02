import PhaseTools.Registry

namespace PhaseTools

inductive NoCallReason where
  | assayNotObservable
  | missingValidatedEnrichment
  | insufficientEvidence
  | ambiguousTopScore
deriving DecidableEq, Repr

inductive CallStatus where
  | called
  | noCall (reason : NoCallReason)
deriving DecidableEq, Repr

structure CallWitness where
  target : Target
  assay : Assay
  validatedEnrichment : Bool
  backend : Backend
  status : CallStatus
  callCount : Nat
deriving Repr

def assayDeclarationValid (witness : CallWitness) : Prop :=
  match witness.assay with
  | .wgs => witness.validatedEnrichment = false
  | .wes => True

def assayFailure
    (target : Target)
    (assay : Assay)
    (validatedEnrichment : Bool) : Option NoCallReason :=
  match assay with
  | .wgs => none
  | .wes =>
      match wesSupport target with
      | .standardCapture => none
      | .validatedTargetEnrichment =>
          match validatedEnrichment with
          | true => none
          | false => some .missingValidatedEnrichment
      | .unsupported => some .assayNotObservable

def evidenceNoCallAllowed (target : Target) (reason : NoCallReason) : Bool :=
  match implementationFor target, reason with
  | .runnable, .insufficientEvidence => true
  | .solverKernel, .insufficientEvidence => true
  | .solverKernel, .ambiguousTopScore => true
  | _, _ => false

def callInvariant (witness : CallWitness) : Prop :=
  assayDeclarationValid witness ∧
    witness.backend = backendFor witness.target ∧
    match witness.status with
    | .called =>
        implemented witness.target = true ∧
          assayFailure witness.target witness.assay witness.validatedEnrichment = none ∧
          0 < witness.callCount
    | .noCall reason =>
        witness.callCount = 0 ∧
          match assayFailure witness.target witness.assay witness.validatedEnrichment with
          | some expected => reason = expected
          | none => evidenceNoCallAllowed witness.target reason = true

instance callInvariantDecidable (witness : CallWitness) :
    Decidable (callInvariant witness) := by
  unfold callInvariant
  cases hAssay : witness.assay <;>
    cases hStatus : witness.status <;>
      cases hFailure : assayFailure witness.target witness.assay
        witness.validatedEnrichment <;>
        simp [assayDeclarationValid, hAssay, hStatus, hFailure] <;>
        infer_instance

def verifyCall (witness : CallWitness) : Bool :=
  decide (callInvariant witness)

theorem verifyCall_sound (witness : CallWitness)
    (verified : verifyCall witness = true) :
    callInvariant witness := by
  have decided : decide (callInvariant witness) = true := by
    simpa [verifyCall] using verified
  exact of_decide_eq_true decided

theorem assayFailure_none_iff_observable
    (target : Target)
    (assay : Assay)
    (validatedEnrichment : Bool) :
    assayFailure target assay validatedEnrichment = none ↔
      observable target assay validatedEnrichment = true := by
  cases target <;> cases assay <;> cases validatedEnrichment <;> rfl

theorem verified_backend_matches_target (witness : CallWitness)
    (verified : verifyCall witness = true) :
    witness.backend = backendFor witness.target := by
  exact (verifyCall_sound witness verified).2.1

theorem verified_called_is_implemented (witness : CallWitness)
    (status : witness.status = .called)
    (verified : verifyCall witness = true) :
    implemented witness.target = true := by
  have invariant := verifyCall_sound witness verified
  rw [status] at invariant
  exact invariant.2.2.1

theorem verified_called_is_observable (witness : CallWitness)
    (status : witness.status = .called)
    (verified : verifyCall witness = true) :
    observable witness.target witness.assay witness.validatedEnrichment = true := by
  have invariant := verifyCall_sound witness verified
  rw [status] at invariant
  exact (assayFailure_none_iff_observable
    witness.target witness.assay witness.validatedEnrichment).mp invariant.2.2.2.1

theorem verified_noCall_has_zero_calls (witness : CallWitness)
    (reason : NoCallReason)
    (status : witness.status = .noCall reason)
    (verified : verifyCall witness = true) :
    witness.callCount = 0 := by
  have invariant := verifyCall_sound witness verified
  rw [status] at invariant
  exact invariant.2.2.1

structure SelectionWitness where
  candidateCount : Nat
  winnerIndex : Nat
  winnerPenalty : Nat
  runnerUpPenalty : Nat
  requiredMargin : Nat
deriving Repr

def selectionInvariant (witness : SelectionWitness) : Prop :=
  0 < witness.candidateCount ∧
    witness.winnerIndex < witness.candidateCount ∧
    witness.winnerPenalty + witness.requiredMargin ≤ witness.runnerUpPenalty

instance selectionInvariantDecidable (witness : SelectionWitness) :
    Decidable (selectionInvariant witness) := by
  unfold selectionInvariant
  infer_instance

def verifySelection (witness : SelectionWitness) : Bool :=
  decide (selectionInvariant witness)

theorem verifySelection_sound (witness : SelectionWitness)
    (verified : verifySelection witness = true) :
    selectionInvariant witness := by
  have decided : decide (selectionInvariant witness) = true := by
    simpa [verifySelection] using verified
  exact of_decide_eq_true decided

theorem verified_winner_in_bounds (witness : SelectionWitness)
    (verified : verifySelection witness = true) :
    witness.winnerIndex < witness.candidateCount := by
  exact (verifySelection_sound witness verified).2.1

theorem verified_margin (witness : SelectionWitness)
    (verified : verifySelection witness = true) :
    witness.winnerPenalty + witness.requiredMargin ≤ witness.runnerUpPenalty := by
  exact (verifySelection_sound witness verified).2.2

end PhaseTools
