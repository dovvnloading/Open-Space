import Lean

/-! Exact models of selected Open-Space mathematical operations.
    Integer time units and scaled natural numbers are intentional:
    no floating-point or astronomical modeling assumptions are hidden here. -/
namespace OpenSpace

def ttMinusUtcUs (taiMinusUtcSeconds : Int) : Int :=
  32184000 + 1000000 * taiMinusUtcSeconds

theorem nicerClockOffset : ttMinusUtcUs 37 = 69184000 := by decide

theorem nicer2018PrintedLabels :
    (43163609441 : Int) - 43094425441 = ttMinusUtcUs 37 := by decide

theorem nicer2022PrintedLabels :
    (86416457212 : Int) - 86347273212 = ttMinusUtcUs 37 := by decide

/-- A common offset never changes a duration or inter-event spacing. -/
theorem commonShiftPreservesDuration (start stop offset : Int) :
    (stop + offset) - (start + offset) = stop - start := by omega

/-- Membership in a closed time window is invariant under a common shift. -/
theorem commonShiftPreservesWindow (start stop t offset : Int) :
    (start ≤ t ∧ t ≤ stop) ↔
      (start + offset ≤ t + offset ∧ t + offset ≤ stop + offset) := by omega

/-- Exact length of an interval intersection (zero if disjoint). -/
def overlapMs (a b c d : Int) : Int :=
  let start := if a ≤ c then c else a
  let stop := if b ≤ d then b else d
  if start ≤ stop then stop - start else 0

theorem orderedOverlap (a b c d : Int)
    (hac : a ≤ c) (hcb : c ≤ b) (hbd : b ≤ d) :
    overlapMs a b c d = b - c := by
  simp [overlapMs, hac, hcb, hbd]

theorem orderedCoveragePartition (a b c d : Int)
    (hac : a ≤ c) (hcb : c ≤ b) (hbd : b ≤ d) :
    (c - a) + overlapMs a b c d + (d - b) = d - a := by
  rw [orderedOverlap a b c d hac hcb hbd]
  omega

/-- In angular units, increasing the modeled limb radius cannot make an
    already covered direction clear, at fixed separation and source offset.
    This is conditional on the supplied limb radii, not a proof about arcsin. -/
theorem largerLimbPreservesCoverage (sep limb₁ limb₂ offset : Int)
    (h : limb₁ ≤ limb₂) (covered : sep - limb₁ - offset ≤ 0) :
    sep - limb₂ - offset ≤ 0 := by omega

/-- Accounting used by open_space/calibrate.py; does not assert recovery
    probabilities or independence of injected trials. -/
theorem calibrationAccounting (attempted recovered missed confounded : Nat)
    (h : attempted = recovered + missed + confounded) :
    attempted - confounded = recovered + missed ∧
      recovered ≤ attempted - confounded := by omega

theorem noAssessableIff (attempted recovered missed confounded : Nat)
    (h : attempted = recovered + missed + confounded) :
    attempted - confounded = 0 ↔ recovered = 0 ∧ missed = 0 := by omega

/-- NICER paired event counts; a missing count is ineligible, not zero. -/
structure Window where
  beforeCount : Option Nat
  afterCount : Option Nat
  otherSource : Bool

structure Summary where
  total : Nat
  eligible : Nat
  wins : Nat
  ties : Nat
  losses : Nat
  beforeEvents : Nat
  afterEvents : Nat
  otherEligible : Nat
  otherWins : Nat
  otherTies : Nat
  otherLosses : Nat
  otherBefore : Nat
  otherAfter : Nat
  deriving DecidableEq, Repr

def zeroSummary : Summary := ⟨0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0⟩

def addSummary (a b : Summary) : Summary :=
  ⟨a.total + b.total, a.eligible + b.eligible, a.wins + b.wins,
   a.ties + b.ties, a.losses + b.losses, a.beforeEvents + b.beforeEvents,
   a.afterEvents + b.afterEvents, a.otherEligible + b.otherEligible,
   a.otherWins + b.otherWins, a.otherTies + b.otherTies,
   a.otherLosses + b.otherLosses, a.otherBefore + b.otherBefore,
   a.otherAfter + b.otherAfter⟩

def classifyWindow (w : Window) : Summary :=
  match w.beforeCount, w.afterCount with
  | some a, some b =>
    let win : Nat := if a < b then 1 else 0
    let tie : Nat := if a = b then 1 else 0
    let loss : Nat := if b < a then 1 else 0
    ⟨1, 1, win, tie, loss, a, b,
     if w.otherSource then 1 else 0,
     if w.otherSource then win else 0,
     if w.otherSource then tie else 0,
     if w.otherSource then loss else 0,
     if w.otherSource then a else 0,
     if w.otherSource then b else 0⟩
  | _, _ => ⟨1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0⟩

def summarize : List Window → Summary
  | [] => zeroSummary
  | w :: ws => addSummary (classifyWindow w) (summarize ws)

/-- Every eligible pair has exactly one outcome. -/
theorem pairedOutcomePartition (a b : Nat) :
    (if a < b then (1 : Nat) else 0) +
    (if a = b then 1 else 0) +
    (if b < a then 1 else 0) = 1 := by
  by_cases h : a < b
  · have hne : a ≠ b := by omega
    have hn : ¬ b < a := by omega
    simp [h, hne, hn]
  · by_cases he : a = b
    · have hn : ¬ b < a := by omega
      simp [h, he, hn]
    · have hl : b < a := by omega
      simp [h, he, hl]

/-- OEC audit: q is a/R_star scaled by 100; r is R_star/R_sun
    scaled by 100. The adopted nominal solar radius is 465046726/10^11 AU.
    Both axis numerators use the same denominator 10^15 AU. -/
def nominalAxisNumerator (q : Nat) : Nat := q * 465046726 * 100
def correctedAxisNumerator (q r : Nat) : Nat := q * 465046726 * r

theorem stellarRadiusScaling (q r : Nat) :
    correctedAxisNumerator q r * 100 = nominalAxisNumerator q * r := by
  unfold correctedAxisNumerator nominalAxisNumerator
  ac_rfl

/-- Exact, division-free checks for the six reported central-value axes.
    o is catalog AU scaled by 1000, rounded is corrected AU scaled by 10^5.
    The 17–23 bound uses F_correct/F_catalog = (a_catalog/a_correct)^2,
    conditional on the same luminosity and the inverse-square model. -/
def axisChecks (q r o rounded : Nat) : Prop :=
  let nominal := nominalAxisNumerator q
  let correct := correctedAxisNumerator q r
  let reported := o * 1000000000000
  let target := rounded * 10000000000
  0 < q ∧ 0 < r ∧
  197 * nominal ≤ 200 * reported ∧ 200 * reported ≤ 203 * nominal ∧
  2 * correct ≤ 2 * target + 10000000000 ∧
  2 * target ≤ 2 * correct + 10000000000 ∧
  4 * correct < reported ∧ reported < 5 * correct ∧
  17 * (correct * correct) < reported * reported ∧
  reported * reported < 23 * (correct * correct)

#print axioms nicerClockOffset
#print axioms nicer2018PrintedLabels
#print axioms nicer2022PrintedLabels
#print axioms commonShiftPreservesDuration
#print axioms commonShiftPreservesWindow
#print axioms orderedOverlap
#print axioms orderedCoveragePartition
#print axioms largerLimbPreservesCoverage
#print axioms calibrationAccounting
#print axioms noAssessableIff
#print axioms pairedOutcomePartition
#print axioms stellarRadiusScaling

end OpenSpace
