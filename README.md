# pg-phenotype against outside oracles and published values

Runs both bindings of [pg-phenotype](https://github.com/rwaples/pg-phenotype)
(the Python wheel and the `pgphenotype` R package) on fixed inputs, and
compares each result with an outside reference implementation (an *oracle*)
or a published value. The latest results are in
[`results/results.md`](results/results.md); every row that is not PASS is
explained in [`differences.md`](differences.md).

## Run it

With [pixi](https://pixi.sh):

```bash
pixi run harness                          # every case (installs the R package first)
pixi run harness --only tet_balanced,krebs
pixi run verify-inputs                    # regenerate inputs/ and compare sha256 with cases.yaml
pixi run lint
```

The first `harness` runs `r-setup`, which downloads the pinned `pgphenotype`
source tarball (see [`pgphenotype.pin`](pgphenotype.pin)), refuses it unless
its sha256 matches, and compiles it into `.r-lib/` (about 45 s; the tarball
builds its vendored Rust crates, so the environment carries `rust`). Later
runs skip the install. A full run takes about 15 s. The script exits 1 if
any row is FAIL or ERROR. There is no CI.

### Checking a local pg-phenotype checkout

Both bindings are pinned to the v0.5.0 release pair: PyPI
`pg-phenotype==0.5.0` in `pixi.toml` and the GitHub release asset
`pgphenotype_0.5.0.tar.gz`. To check a candidate instead, without committing
the change:

1. In `pixi.toml`, replace the `pg-phenotype` entry with
   `pg-phenotype = { path = "../pg-phenotype", editable = true }` and run
   `pixi install`.
2. Build the R tarball in the checkout (`pixi run -e r tools/r_build_tarball.sh`)
   and run `PGP_VALIDATE_R_TARBALL=/path/to/pgphenotype_<ver>.tar.gz pixi run harness`.

The preflight then requires the two bindings to link the same
pedigree-graph-core revision instead of the pinned version. Without the
override it refuses to run unless both report 0.5.0.

## What it does

For each case in [`cases.yaml`](cases.yaml), `run_oracles.py`:

1. **Builds the input.** Contingency tables are literal in `cases.yaml`.
   Pedigrees and survival data are committed in `inputs/*.tsv.gz`, written by
   the seeded [`make_inputs.py`](make_inputs.py), and the harness refuses a
   file whose sha256 differs from `cases.yaml`.
2. **Runs pg-phenotype twice**: the Python binding in-process and the R
   binding through [`pg_r.R`](pg_r.R). The two must agree to 1e-12 on every
   compared value.
3. **Runs the oracles.** Each is an R script in [`oracles/`](oracles) with
   one interface, `Rscript <script> <input.json> <output.json>`, writing
   doubles with 17 significant digits.
4. **Compares** each check at its tolerance, after any declared conversion
   ([`conversions.py`](conversions.py)), and writes `results/results.md` and
   `results/results.tsv`. Each row records the versions of R, the oracle
   package and both bindings.

## Cases

| Family | pg-phenotype | Oracle or source |
|---|---|---|
| Tetrachoric | `correlation.tetrachoric(table=)`, `correlation_tetrachoric()` | `polycor::polychor` (two-step rho, joint-ML SE, fixed-threshold SE) |
| Polychoric | `assortative.mate_correlation` on one-pair Mate Networks | `polycor::polychor`, `psych::polychoric(correct = 0, global = FALSE)` |
| Polyserial, biserial | the same | `polycor::polyserial(ML = FALSE, std.err = TRUE)` |
| Parent-offspring cluster SE | `regression.parent_offspring`, `regression_parent_offspring()` | `lm`, `sandwich::vcovCL(type = "HC0", cadjust = TRUE)` |
| PA-FGRS single-relative values | `pafgrs.score_univariate` on an MZ pair at h² = 1 | closed form (scipy); Dybdahl Krebs et al. (2024), supp. Note S3 |
| Aalen-Johansen CIP | `Cip.from_aalen_johansen`, `pafgrs_cip_aalen_johansen()` | `survival::survfit`, `cmprsk::cuminc` |
| Relative-Pair Correlation SE | not checked (GAP G1) | |

The PA-FGRS source is Dybdahl Krebs M, Georgii Hellberg K-L, Lundberg M, et
al. (2024). Genetic liability estimated from large-scale family data
improves genetic prediction, risk score profiling, and gene mapping for
major depression. *Am J Hum Genet* 111(11):2494-2509.
doi:[10.1016/j.ajhg.2024.09.009](https://doi.org/10.1016/j.ajhg.2024.09.009).

## Oracles

Every oracle is installed from conda-forge and called; none of their code is
copied here. Versions are the ones in `pixi.lock`, licences as each
installed package's `DESCRIPTION` states them.

| Package | Version | Licence | Used for |
|---|---|---|---|
| R (`r-base`) | 4.5.3 | GPL-2.0-or-later (conda metadata) | |
| polycor | 0.8.2 | GPL (>= 2) | tetrachoric, polychoric, polyserial |
| psych | 2.6.9 | GPL (>= 2) | polychoric |
| sandwich | 3.1.3 | GPL-2 \| GPL-3 | cluster-robust SE |
| survival | 3.8.12 | LGPL (>= 2) | Aalen-Johansen |
| cmprsk | 2.2.12 | GPL (>= 2) | cumulative incidence |
| jsonlite | 2.0.0 | MIT | the JSON interface |
| scipy | 1.18.1 | BSD-3-Clause (conda metadata) | closed-form PA-FGRS rows |

## Layout

```
cases.yaml        cases and their checks (schema in its header)
run_oracles.py    the driver: preflight, inputs, both bindings, oracles, report
pg_r.R            the R binding on each case
oracles/*.R       one script per oracle package; _io.R holds the JSON interface
conversions.py    the closed list of declared conversions (C1, C2; C3 reserved)
make_inputs.py    the seeded generator of inputs/*.tsv.gz
differences.md    every DIFF, GAP and INFO row explained
docs/adr/         decisions; GLOSSARY.md the terms
```

## Licence

MIT; see [LICENSE](LICENSE).
