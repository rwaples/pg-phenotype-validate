# Where pg-phenotype's results differ from the references

Notes on the [harness results](results/results.md) at pg-phenotype 0.5.0
(both bindings, pedigree-graph-core `0a79c4e`). Every row whose reference is
the same estimator under the same convention passes: the tetrachoric rho
and SE, the polychoric rho, the parent-offspring slope, OLS SE and Mate
Network sandwich SE, the PA-FGRS closed forms and the Aalen-Johansen CIF.
The rows below are a known difference (D), a check not made (G), or a
different estimator shown for information (I). Conversions C1 and C2 are
listed in [`conversions.py`](conversions.py).

## D1. Tetrachoric SE on a table with an empty cell

`tet_empty_cell`, table `[[45, 0], [15, 20]]`. pg sets `boundary` and
withholds the SE (`se_unavailable_reason = boundary`;
pg-phenotype `crates/core/src/correlation.rs:149-151`). polycor's joint-ML
SE there is 23.36: its numerical Hessian is taken on a likelihood that is
flat towards rho = 1, so the value has no meaning as an SE. The two rhos
agree within polycor's tolerance (0.99764 against 0.99761).

## D2. The published PA-FGRS case value, 2.07

Dybdahl Krebs et al. (2024), supp. p.10 (Note S3), print 2.07 for a case at
K = 0.05 at every onset age. The closed form E[L | L > T] = φ(T)/K is
2.0627128 (T = 1.6448536), which pg reproduces to 4e-16 (the closed-form
row), and 2.0627 does not round to 2.07 (|diff| 0.0073, tolerance 0.005,
half a unit in the last printed place).

We take 2.07 to be a slip in the table, for these reasons:

- Every other cell of the table reproduces at its printed precision:
  - the simple-liability control at k = 0.05 (−0.1086 against −0.11);
  - the simple-liability case and control at k = 0.01 (2.6652 and −0.0269
    against 2.67 and −0.03);
  - the mixture controls at k = 0.05, observed to 60, 30 and 1 (−0.10856,
    −0.05289 and −0.00172 against −0.11, −0.05 and −0.002), and at
    k = 0.01 (−0.02084, −0.01037 and −0.000344 against −0.02, −0.01 and
    −0.0003), from −w φ(T)/(1 − wK) with w = CIP(age)/K;
  - the ADT cases at onset 60, 30 and 1, at k = 0.05 (1.6449, 1.9600,
    3.1440 against 1.64, 1.96, 3.14) and k = 0.01 (2.3263, 2.5758, 3.5879
    against 2.33, 2.58, 3.59). These are computed as Φ⁻¹(1 − CIP(age)), a
    reading of the ADT model inferred from the match and not checked
    against the paper's definition.
- The off value is a single number. All nine cells that should hold φ(T)/K
  print 2.07: the simple model at k = 0.05 and the mixture model at both
  k, three onset ages each. The mixture case does not depend on k.
- Nothing natural produces 2.07:
  - PAFGRS's `FGRS()` integrates numerically,
    `integrate(qnorm, 1 - prev, 1) / prev` (PAFGRS `R/est_liab.R:173`,
    commit `d272e32`). At K = 0.05 that gives 2.062713, with an
    absolute error of 4e-7.
  - Rounding T to 1.64 gives 2.0792, which prints as 2.08. Printing 2.07
    needs T in [1.6412, 1.6442].
  - Printing 2.07 at T = Φ⁻¹(1 − K) needs K in [0.0485, 0.0497], while the
    controls fit K = 0.05.
  - PAFGRS's scoring uses the same closed-form truncated-normal mean
    (`R/est_liab.R:63`), so it gives 2.0627 too.

The note does not say how the table was computed, so a route to 2.07 not
listed here (a simulation, say) cannot be ruled out. Only the authors could
confirm. The rows stay DIFF, since the published value is what they test.
Figures recomputed in R 4.5.3 on 2026-10-21.

## D3. Polyserial and biserial rho: the standardisation divisor

`polyserial` and `biserial`, 500 pairs each.

| Case | polycor | pg | gap |
|---|---|---|---|
| polyserial | 0.3852798 | 0.3849702 | 3.1e-4 |
| biserial | 0.2850178 | 0.2847595 | 2.6e-4 |

pg standardises the continuous side with the 1/N variance
(pg-phenotype `docs/assortative-mating.md`, Estimators); polycor uses
`scale()`, the n − 1 one. [`docs/c3_divisor.R`](docs/c3_divisor.R)
maximises the two-step likelihood (Olsson, Drasgow & Dorans 1982, eqs.
19–20), with polycor's first-step thresholds, by `optimize()`. The
likelihood is flat in rho at about 1e-8 in double precision, so these fits
are good to about 1e-8, not to `optimize()`'s tol of 1e-14:

- with `scale(x)` it gives 0.385279642 and 0.285017785, so polycor's BFGS
  is within 2e-7 of its optimum and the gap is not its precision;
- with the 1/N standardisation it gives 0.384970228 and 0.284759538,
  within 1e-8 of pg, which is as close as these fits can show.

The divisor explains the whole gap, but not by a fixed factor: the ratio of
the two fits is 0.9991969 in one case and 0.9990939 in the other, and
sqrt((n − 1)/n) is 0.9989995. Rescaling the probit slope ρ/√(1 − ρ²) by
sqrt((n − 1)/n) instead predicts pg's values only to 1.9e-5 and 3.8e-6,
because the thresholds stay at the marginal values while z is rescaled. So
C3 stays reserved and these rows are DIFF. The gap is O(1/n), against pg SEs
of 0.039 and 0.062 (the `SE vs polycor two-step SE` rows of
[`results/results.md`](results/results.md)).

pg's divisor is the definition's. Olsson, Drasgow & Dorans (1982; full
reference below) define
the two-step estimator by estimating the mean and variance of X with
F₃ = σ² − (1/N) Σ (xᵢ − X̄)² = 0 (eq. 36, p.344), the 1/N variance, before
maximising the likelihood (eq. 20, p.342) in ρ alone. polycor's `scale()`
is a choice of its own. pg-phenotype standardises with 1/N
(pg-phenotype `docs/assortative-mating.md`, Estimators), so D3 is accepted
([`status.md`](status.md)). It reopens if pg-phenotype changes that
divisor. Equations and pages checked against the paper's PDF on
2026-10-21.

Olsson U, Drasgow F, Dorans NJ (1982). The polyserial correlation
coefficient. *Psychometrika* 47(3):337–347. doi:10.1007/BF02294164.

## G1. Relative-Pair Correlation SE

pg-phenotype's Relative-Pair Correlation reports no per-pair or per-entry
influence (pairs are never stored; pg-phenotype ADR 0013), so its two-way
cluster SE cannot be rebuilt by an outside sandwich. The row waits on a
pg-phenotype issue, [rwaples/pg-phenotype#25](https://github.com/rwaples/pg-phenotype/issues/25),
asking for opt-in per-entry influence output (member ids,
MRCA key, influence).

## G2. The PA-FGRS mixture column at k = 0.01

Note S3's mixture column at k = 0.01 takes the threshold from the true
K = 0.05 and the weights from an observed CIP that tops out at 0.01. pg's
`Cip` takes K from its last value (`pg_phenotype/pafgrs.py`, `Cip`), so it
cannot express a threshold apart from the table. The ADT column is not a pg
method and is out of scope.

## G3. cuminc with late entry

`cmprsk::cuminc(ftime, fstatus, group, strata, rho = 0, cencode = 0,
subset, na.action = na.omit)` has no entry argument, so `aj_late_entry` is
checked against `survfit` only.

## INFO rows: different estimators

Shown with their measured gap (relative after C2 for the mate-correlation
rows); not counted.

- **Tetrachoric, fixed-threshold SE.** `polychor(ML = FALSE, std.err = TRUE)`
  holds the thresholds fixed. On every table it is below pg's SE, as it
  should be: by 1.2e-5 (balanced, n = 1000) to 3.4e-3 (a cell of 1, n = 80),
  absolute.
- **Polychoric SE.** pg's two-step sandwich against polycor's joint-ML SE:
  3.5% (`poly_3x4`) and 0.15% (`poly_4x4`) relative; against the
  fixed-threshold SE, 5.1% and 1.0%. Codex's review reported 5% and 10% on
  two other tables, so the size depends on the table.
- **Polyserial and biserial SE.** polycor's two-step SE is the inverse
  Hessian with thresholds fixed, not a sandwich: 0.9% and 7.6% relative.

## Notes on the plan

- jsonlite's `digits = NA` writes 15 significant digits, which does not
  round-trip a double; the scripts write `digits = I(17)`
  (`oracles/_io.R`).
- polycor's two-step tetrachoric rho on the balanced table is 8.57e-7 from
  pg's, as Codex measured (8.6e-7).
