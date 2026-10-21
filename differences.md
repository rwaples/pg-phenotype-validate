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
half a unit in the last printed place). The table's other cells at the same
K agree within that tolerance (−0.11, −0.05, −0.002 against −0.1086,
−0.0529, −0.0017), and its k = 0.01 simple-liability cells too (2.67 and
−0.03 against 2.6652 and −0.0269). What produced 2.07 is not established.

## D3. Polyserial and biserial rho: the standardisation divisor

`polyserial` and `biserial`, 500 pairs each.

| Case | polycor | pg | gap |
|---|---|---|---|
| polyserial | 0.3852798 | 0.3849702 | 3.1e-4 |
| biserial | 0.2850178 | 0.2847595 | 2.6e-4 |

pg standardises the continuous side with the 1/N variance
(pg-phenotype `docs/assortative-mating.md`, Estimators); polycor uses
`scale()`, the n − 1 one. [`docs/c3_divisor.R`](docs/c3_divisor.R) refits
polycor's own two-step likelihood with `optimize(tol = 1e-14)`:

- with `scale(x)` it gives 0.385279654 and 0.285017776, so polycor's BFGS
  is within 2e-7 of its optimum and the gap is not its precision;
- with the 1/N standardisation it gives 0.384970222 and 0.284759526,
  within 1.3e-8 of pg.

The divisor explains the whole gap, but not by a fixed factor: the ratio of
the two fits is 0.9991969 in one case and 0.9990939 in the other, and
sqrt((n − 1)/n) is 0.9989995. So C3 stays reserved and these rows are DIFF.
The gap is O(1/n).

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
