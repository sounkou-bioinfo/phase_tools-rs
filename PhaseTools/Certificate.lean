namespace PhaseTools

inductive CallStatus where
  | called
  | noCall
deriving DecidableEq, Repr

structure CallWitness where
  status : CallStatus
  callCount : Nat
deriving Repr

def callInvariant (witness : CallWitness) : Prop :=
  match witness.status with
  | .called => 0 < witness.callCount
  | .noCall => witness.callCount = 0

instance callInvariantDecidable (witness : CallWitness) :
    Decidable (callInvariant witness) := by
  unfold callInvariant
  cases witness.status <;> infer_instance

def verifyCall (witness : CallWitness) : Bool :=
  decide (callInvariant witness)

theorem verifyCall_sound (witness : CallWitness)
    (verified : verifyCall witness = true) :
    callInvariant witness := by
  have decided : decide (callInvariant witness) = true := by
    simpa [verifyCall] using verified
  exact of_decide_eq_true decided

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
