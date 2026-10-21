# ADR 0001: Oracles are called, and convention gaps are declared

**Status:** accepted
**Date:** 2026-10-21

## Context

pg-phenotype's own tests check it against golden values it generated, some
from the same R packages. A separate public repo can check both bindings
against outside implementations at pinned versions, as pedsum-validate does
for pedsum. Two questions decide whether its rows can be trusted: where the
reference numbers come from, and what happens when an oracle's convention
differs from pg-phenotype's.

The oracles (polycor, psych, sandwich, cmprsk) are GPL; this repo is MIT.
Several known gaps are pure conventions: polycor's SE lacks pg's
n/(n-1) (C1) or G/(G-1) (C2) factor. Others are different estimators: on an
r x s table, polycor's joint-ML SE is not pg's two-step sandwich SE.

## Decision

1. **Oracles are installed and called, never vendored.** Each is a
   conda-forge package at the version `pixi.lock` fixes, run through an R
   script with one JSON interface. Oracle values live only here;
   pg-phenotype's `tests/golden/` is unchanged.
2. **A tolerance is the oracle's own precision**, stated per check. polycor
   fits rho with `optimise()` at its default tolerance, so its rho rows run
   at 2e-4; exact code paths (OLS, the sandwich, Aalen-Johansen) run at
   1e-10 to 1e-12.
3. **A convention gap is a Conversion from a closed list** (conversions.py),
   applied to the oracle's value and shown in the row. An unknown id is an
   ERROR, and a gap no Conversion explains is a FAIL.
4. **A different estimator is INFO**, reported with its measured gap and
   never counted as PASS or FAIL.
5. **A known difference that is not a factor is DIFF**, with an entry in
   differences.md. A DIFF that matches becomes a FAIL so the entry is
   revisited.
6. **Both bindings are checked at one release pair** (PyPI wheel and R
   tarball), and the preflight refuses to run when they report another
   version, or, under the local override, different core revisions.

## Consequences

- A FAIL means pg-phenotype moved, or an oracle did. Rows record every
  version, so a relock that moves an oracle shows which.
- A Conversion needs an exact factor. C3 (the 1/N vs n-1 standardisation of
  polyserial) has none, so those rows are DIFF (differences.md D3).
- The harness depends on GPL packages at run time only; nothing GPL is
  distributed.
