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

def assayDeclarationValid (witness : CallWitness) : Bool :=
  match witness.assay with
  | .wgs => !witness.validatedEnrichment
  | .wes => true

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

def noAssayFailure
    (target : Target)
    (assay : Assay)
    (validatedEnrichment : Bool) : Bool :=
  match assayFailure target assay validatedEnrichment with
  | none => true
  | some _ => false

def evidenceNoCallAllowed (target : Target) (reason : NoCallReason) : Bool :=
  match implementationFor target, reason with
  | .runnable, .insufficientEvidence => true
  | .solverKernel, .insufficientEvidence => true
  | .solverKernel, .ambiguousTopScore => true
  | _, _ => false

def callStatusValid (witness : CallWitness) : Bool :=
  match witness.status with
  | .called =>
      implemented witness.target &&
        (noAssayFailure witness.target witness.assay witness.validatedEnrichment &&
          decide (0 < witness.callCount))
  | .noCall reason =>
      decide (witness.callCount = 0) &&
        (match assayFailure witness.target witness.assay witness.validatedEnrichment with
        | some expected => decide (reason = expected)
        | none => evidenceNoCallAllowed witness.target reason)

def callInvariant (witness : CallWitness) : Prop :=
  assayDeclarationValid witness = true ∧
    witness.backend = backendFor witness.target ∧
    callStatusValid witness = true

def verifyCall (witness : CallWitness) : Bool :=
  assayDeclarationValid witness &&
    (decide (witness.backend = backendFor witness.target) && callStatusValid witness)

theorem verifyCall_sound (witness : CallWitness)
    (verified : verifyCall witness = true) :
    callInvariant witness := by
  simp only [verifyCall, Bool.and_eq_true] at verified
  exact ⟨verified.1, of_decide_eq_true verified.2.1, verified.2.2⟩

theorem noAssayFailure_iff_observable
    (target : Target)
    (assay : Assay)
    (validatedEnrichment : Bool) :
    noAssayFailure target assay validatedEnrichment = true ↔
      observable target assay validatedEnrichment = true := by
  cases target <;> cases assay <;> cases validatedEnrichment <;>
    simp [noAssayFailure, assayFailure, observable, wesSupport]

theorem verified_assay_declaration_valid (witness : CallWitness)
    (verified : verifyCall witness = true) :
    assayDeclarationValid witness = true := by
  exact (verifyCall_sound witness verified).1

theorem verified_backend_matches_target (witness : CallWitness)
    (verified : verifyCall witness = true) :
    witness.backend = backendFor witness.target := by
  exact (verifyCall_sound witness verified).2.1

theorem verified_called_is_implemented (witness : CallWitness)
    (status : witness.status = .called)
    (verified : verifyCall witness = true) :
    implemented witness.target = true := by
  have statusValid := (verifyCall_sound witness verified).2.2
  unfold callStatusValid at statusValid
  rw [status] at statusValid
  simp only [Bool.and_eq_true] at statusValid
  exact statusValid.1

theorem verified_called_is_observable (witness : CallWitness)
    (status : witness.status = .called)
    (verified : verifyCall witness = true) :
    observable witness.target witness.assay witness.validatedEnrichment = true := by
  have statusValid := (verifyCall_sound witness verified).2.2
  unfold callStatusValid at statusValid
  rw [status] at statusValid
  simp only [Bool.and_eq_true] at statusValid
  exact (noAssayFailure_iff_observable
    witness.target witness.assay witness.validatedEnrichment).mp statusValid.2.1

theorem verified_called_has_positive_calls (witness : CallWitness)
    (status : witness.status = .called)
    (verified : verifyCall witness = true) :
    0 < witness.callCount := by
  have statusValid := (verifyCall_sound witness verified).2.2
  unfold callStatusValid at statusValid
  rw [status] at statusValid
  simp only [Bool.and_eq_true] at statusValid
  exact of_decide_eq_true statusValid.2.2

theorem verified_noCall_has_zero_calls (witness : CallWitness)
    (reason : NoCallReason)
    (status : witness.status = .noCall reason)
    (verified : verifyCall witness = true) :
    witness.callCount = 0 := by
  have statusValid := (verifyCall_sound witness verified).2.2
  unfold callStatusValid at statusValid
  rw [status] at statusValid
  simp only [Bool.and_eq_true] at statusValid
  exact of_decide_eq_true statusValid.1

theorem verified_noCall_reason_valid (witness : CallWitness)
    (reason : NoCallReason)
    (status : witness.status = .noCall reason)
    (verified : verifyCall witness = true) :
    (match assayFailure witness.target witness.assay witness.validatedEnrichment with
    | some expected => decide (reason = expected)
    | none => evidenceNoCallAllowed witness.target reason) = true := by
  have statusValid := (verifyCall_sound witness verified).2.2
  unfold callStatusValid at statusValid
  rw [status] at statusValid
  simp only [Bool.and_eq_true] at statusValid
  exact statusValid.2

theorem verified_noCall_exact_assay_reason (witness : CallWitness)
    (reason expected : NoCallReason)
    (status : witness.status = .noCall reason)
    (failure : assayFailure witness.target witness.assay witness.validatedEnrichment =
      some expected)
    (verified : verifyCall witness = true) :
    reason = expected := by
  have reasonValid := verified_noCall_reason_valid witness reason status verified
  rw [failure] at reasonValid
  exact of_decide_eq_true reasonValid

theorem verified_noCall_evidence_reason_allowed (witness : CallWitness)
    (reason : NoCallReason)
    (status : witness.status = .noCall reason)
    (noFailure : assayFailure witness.target witness.assay witness.validatedEnrichment = none)
    (verified : verifyCall witness = true) :
    evidenceNoCallAllowed witness.target reason = true := by
  have reasonValid := verified_noCall_reason_valid witness reason status verified
  rw [noFailure] at reasonValid
  exact reasonValid

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
