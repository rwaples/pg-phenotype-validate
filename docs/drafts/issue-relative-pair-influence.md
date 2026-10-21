<!-- Draft for rwaples/pg-phenotype. Title: Opt-in per-entry influence output for the Relative-Pair Correlation SE -->

## What

An opt-in output from `relative_pairs.correlation(..., se=True)` (and
`relative_pairs_correlation()` in R) that returns, for each cell, one record
per pair entry: the two member ids, the pair's MRCA key, and the entry's
influence φ on r̂.

## Why

pg-phenotype-validate checks pg-phenotype's SEs against outside sandwich
implementations. The parent-offspring Mate Network SE is checked that way
against `sandwich::vcovCL` to 1e-10. The Relative-Pair Correlation SE
cannot be: it clusters two pairs when they share a member or an MRCA, by
inclusion-exclusion (ADR 0013, Decision: `meat = Σ_m Φ_m² + (Σ_i S_i² −
Σ_d φ_d²) − (Σ_(i,m) S_(i,m)² − Σ_d φ_d²)`), and pairs are never stored
(ADR 0013:69), so nothing outside the core sees the φ_d or the cluster keys.

With per-entry ids, MRCA keys and influences, an outside script can rebuild
the meat (for example as a two-way cluster sum, `sandwich::vcovCL` with two
cluster variables or a direct inclusion-exclusion over the three sums) and
check the SE. Today pg-phenotype-validate reports this as GAP G1.

## Shape (a suggestion)

- Off by default, so memory stays as it is (pairs are never stored unless asked).
- Per cell: `entries` with `first_id`, `second_id`, `mrca_id` (the row ADR
  0013 picks; a founder for a founder MRCA; the pair itself when there is
  none), and `influence`; for a symmetric pair, both entries or one record
  with the averaged φ, whichever the SE uses.
- The same names in Python and R.
