# pg-phenotype-validate

Checks pg-phenotype's Python and R bindings against outside oracles and
published values. Terms are in GLOSSARY.md, decisions in docs/adr/.

- Everything runs through pixi: `pixi run harness`, `pixi run lint`,
  `pixi run verify-inputs`. `pixi.lock` changes only on purpose.
- Never copy oracle code (GPL); call it from `oracles/*.R` (ADR 0001).
- Tolerances are the oracle's precision; never widen one to make a row
  pass. A new convention gap needs an exact factor in conversions.py, or a
  DIFF entry in differences.md.
- Every non-PASS row has an entry in status.md; update it when a run
  changes a row's result or a follow-up lands.
- After changing a generator, `pixi run make-inputs` and update each
  sha256 in cases.yaml.
- Don't commit the local-checkout override in pixi.toml.
- Before reporting, quote `pixi run lint`, `pixi run harness` with its exit
  code, and the totals line of results/results.md.
