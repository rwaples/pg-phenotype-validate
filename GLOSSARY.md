# pg-phenotype-validate glossary

**Oracle**:
An outside implementation of a method pg-phenotype also implements, installed
and called at a pinned version, never copied (ADR 0001). Its output is the
reference a Check compares with. A published value or a closed form computed
by the driver is a reference too, but not an Oracle.
_Avoid_: golden, ground truth

**Case**:
One input (a table, a pedigree, survival data, or a set of PA-FGRS
scenarios) run through both pg-phenotype bindings and every Oracle its
Checks name. One entry under `cases:` in cases.yaml.

**Check**:
One comparison of one pg-phenotype result with one reference, at a stated
tolerance, giving one row of results. Both bindings must meet it, and agree
with each other to 1e-12.

**Conversion**:
A declared factor that turns an Oracle's value into pg-phenotype's
convention, from the closed list in conversions.py (C1, C2). The row shows
the raw value, the factor and the converted value. A gap that no listed
Conversion explains is a FAIL, not a Conversion. A different estimator is
not a Conversion either; its rows are INFO.
_Avoid_: correction, adjustment

**Status**:
The outcome of a Check: PASS (within tolerance), DIFF (a known difference,
explained under its id in differences.md), GAP (not checked, with an id),
INFO (a different estimator, its gap shown and never counted as PASS or
FAIL), FAIL, ERROR. A FAIL or ERROR makes the harness exit 1. A declared
DIFF that matches is a FAIL, so a stale declaration surfaces.

**Release pair**:
The PyPI wheel and the R source tarball of one pg-phenotype release, the
pin both bindings are checked at (pixi.toml, pgphenotype.pin).
