# Status of every row that is not PASS

One entry per non-PASS row group in [`results/results.md`](results/results.md),
with where it stands and what would close it. Why each one differs is in
[`differences.md`](differences.md), under the same id; this file does not repeat
it. Update the entry when its status changes, and add one when a run produces a
new non-PASS row.

Run: pg-phenotype 0.5.0 (both bindings, core `0a79c4e72f12`), 2026-10-21.
Totals then: DIFF 6, GAP 3, INFO 11, PASS 64.

Statuses:

- **Open**: needs work or a decision; the next step says which.
- **Blocked**: waits on something outside this repo.
- **Accepted**: understood, and no change is planned.
- **Info**: a different estimator, reported for comparison only.

| Id | Rows | Result status | Status | Next step | Updated |
|---|---|---|---|---|---|
| D1 | `tet_empty_cell` SE vs joint-ML SE | DIFF | Accepted | None. pg withholds the SE at a boundary, which is intended; polycor's 23.36 is not a usable SE there. | 2026-10-21 |
| D2 | `krebs` case at onset 60, 30, 1 (K = 0.05) | DIFF ×3 | Open | Find what produced the published 2.07 when the closed form gives 2.0627. One way is to run R PAFGRS (`external/PAFGRS`, commit `d272e32`) on the same one-person setup; if it also gives 2.0627, record the published value as a misprint. | 2026-10-21 |
| D3 | `polyserial`, `biserial` rho | DIFF ×2 | Open | Decide whether pg-phenotype keeps the 1/N standardisation or moves to n − 1 to match polycor. The gap is O(1/n) and is fully explained (`docs/c3_divisor.R`), but there is no fixed factor, so C3 cannot close it. If pg keeps 1/N, mark this Accepted. | 2026-10-21 |
| G1 | `relative_pair_se` | GAP | Blocked | Waits on [rwaples/pg-phenotype#25](https://github.com/rwaples/pg-phenotype/issues/25) (per-entry influence output). Once released, add the vcovCL two-way check. | 2026-10-21 |
| G2 | `krebs` k = 0.01 mixture column | GAP | Accepted | None while pg's `Cip` takes K from its last value. Reopen if pg-phenotype adds a threshold separate from the CIP. | 2026-10-21 |
| G3 | `aj_late_entry` vs cuminc | GAP | Accepted | None. cuminc has no entry argument, and `survfit` covers the late-entry case. | 2026-10-21 |
| I1 | tetrachoric SE vs fixed-threshold SE, 5 tables | INFO ×5 | Info | None. A negative control: pg's SE is above the fixed-threshold SE on every table, as expected. Flag it if that ever flips. | 2026-10-21 |
| I2 | `poly_3x4`, `poly_4x4` SE vs joint-ML and fixed-threshold SE | INFO ×4 | Info | None; no outside oracle computes pg's two-step sandwich on an r × s table. | 2026-10-21 |
| I3 | `polyserial`, `biserial` SE vs polycor two-step SE | INFO ×2 | Info | None; polycor's SE is an inverse Hessian with the thresholds fixed, not a sandwich. | 2026-10-21 |

Open: D2, D3. Blocked: G1.
