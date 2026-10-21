"""Check pg-phenotype's two bindings against outside oracles and published values.

For each case in cases.yaml: build its input, run the Python binding here and
the R binding through pg_r.R, run each oracle's R script, and compare.  Writes
results/results.md and results/results.tsv; exits 1 on any FAIL or ERROR.

    python run_oracles.py               # every case
    python run_oracles.py --only tet_balanced,krebs
    python run_oracles.py --install-r   # install the pinned R tarball into .r-lib/
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
import urllib.request
from dataclasses import dataclass, field
from importlib.metadata import version as dist_version
from pathlib import Path
from typing import Any

import numpy as np
import scipy
import yaml
from scipy.stats import norm

from conversions import CONVERSIONS, RESERVED

ROOT = Path(__file__).resolve().parent
WORK = ROOT / ".work"
R_LIB = ROOT / ".r-lib"
RESULTS = ROOT / "results"
#: Set to a local pgphenotype source tarball to check a checkout (not committed).
OVERRIDE_ENV = "PGP_VALIDATE_R_TARBALL"
#: The two bindings must agree to this, on every compared value.
BINDINGS_TOL = 1e-12
STATUSES = ("PASS", "DIFF", "GAP", "INFO", "FAIL", "ERROR")
FAMILIES = ("tetrachoric", "mate", "parent_offspring", "krebs", "aalen_johansen", "gap")
LINES = ("midparent", "mother", "father")


# ---------------------------------------------------------------- R setup


def read_pin() -> dict[str, str]:
    """pgphenotype.pin: the R tarball's version, url and sha256."""
    return yaml.safe_load((ROOT / "pgphenotype.pin").read_text())


def sha256(path: Path) -> str:
    """The file's sha256, hex."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def installed_r_version() -> str | None:
    """The pgphenotype version in .r-lib/, or None."""
    description = R_LIB / "pgphenotype" / "DESCRIPTION"
    if not description.exists():
        return None
    for line in description.read_text().splitlines():
        if line.startswith("Version:"):
            return line.split(":", 1)[1].strip()
    return None


def install_r() -> None:
    """Install the pinned tarball (or the override) into .r-lib/, once per tarball."""
    pin = read_pin()
    override = os.environ.get(OVERRIDE_ENV)
    stamp = R_LIB / ".installed-sha256"
    if override:
        tarball = Path(override).resolve()
        want = sha256(tarball)
    else:
        tarball = R_LIB / "_download" / f"pgphenotype_{pin['version']}.tar.gz"
        want = pin["sha256"]
    if stamp.exists() and stamp.read_text().strip() == want and installed_r_version():
        print(f"r-setup: pgphenotype {installed_r_version()} already installed from {want[:12]}")
        return
    if not override and not (tarball.exists() and sha256(tarball) == want):
        tarball.parent.mkdir(parents=True, exist_ok=True)
        print(f"r-setup: downloading {pin['url']}")
        with urllib.request.urlopen(pin["url"]) as response, tarball.open("wb") as out:
            shutil.copyfileobj(response, out)
    got = sha256(tarball)
    if got != want:
        sys.exit(f"r-setup: {tarball} has sha256 {got}, pinned {want}; refusing it")
    R_LIB.mkdir(exist_ok=True)
    subprocess.run(["R", "CMD", "INSTALL", f"--library={R_LIB}", str(tarball)], check=True)
    stamp.write_text(want + "\n")
    print(f"r-setup: installed pgphenotype {installed_r_version()} from {want[:12]}")


# ---------------------------------------------------------------- inputs


def read_input(spec: dict[str, str]) -> dict[str, list[Any]]:
    """A committed inputs/*.tsv.gz as columns, refusing a file whose sha256 differs from the pin."""
    path = ROOT / spec["file"]
    got = sha256(path)
    if got != spec["sha256"]:
        sys.exit(f"{spec['file']} has sha256 {got}, cases.yaml pins {spec['sha256']}; refusing it")
    with gzip.open(path, "rt", newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    columns: dict[str, list[Any]] = {name: [] for name in rows[0]}
    for row in rows:
        for name, raw in row.items():
            columns[name].append(_parse(raw))
    return columns


def _parse(raw: str) -> Any:
    if raw == "":
        return None
    try:
        return int(raw)
    except ValueError:
        pass
    try:
        return float(raw)
    except ValueError:
        return raw


def _floats(values: list[Any]) -> np.ndarray:
    return np.array([math.nan if v is None else float(v) for v in values])


# ---------------------------------------------------------------- payloads


@dataclass
class Case:
    """One case of cases.yaml with its built inputs."""

    spec: dict[str, Any]
    payload: dict[str, Any] = field(default_factory=dict)  # the pg input, for both bindings
    oracle_payload: dict[str, Any] = field(default_factory=dict)  # what each oracle script reads
    context: dict[str, Any] = field(default_factory=dict)  # n, G, ... for conversions and checks

    @property
    def id(self) -> str:
        """The case id."""
        return self.spec["id"]

    @property
    def family(self) -> str:
        """The case family."""
        return self.spec["family"]


def build(spec: dict[str, Any]) -> Case:
    """Read a case's input and build what pg and the oracles take."""
    case = Case(spec)
    family = spec["family"]
    if family not in FAMILIES:
        raise ValueError(f"unknown family {family!r}")
    if family == "tetrachoric":
        table = spec["table"]
        case.payload = {"table": table}
        case.oracle_payload = {"table": table}
        case.context = {"n": sum(map(sum, table))}
    elif family == "mate":
        _build_mate(case, read_input(spec["input"]))
    elif family == "parent_offspring":
        _build_po(case, read_input(spec["input"]))
    elif family == "krebs":
        case.payload = {"scenarios": spec["scenarios"]}
    elif family == "aalen_johansen":
        cols = read_input(spec["input"])
        case.payload = {"age": cols["age"], "status": cols["status"], "entry": cols.get("entry")}
        case.oracle_payload = case.payload
    return case


def _build_mate(case: Case, cols: dict[str, list[Any]]) -> None:
    spec = case.spec
    kinds = spec["kinds"]
    traits = [cols[f"t{i}"] for i in range(len(kinds))]
    case.payload = {
        "id": cols["id"],
        "mother": cols["mother"],
        "father": cols["father"],
        "sex": cols["sex"],
        "traits": traits,
        "kinds": kinds,
        "cell": spec["cell"],
    }
    # The oracle's input: (mother's value, father's value) of each Mating Pair.
    row = {rid: i for i, rid in enumerate(cols["id"])}
    i, j = spec["cell"]
    pairs = sorted({(m, f) for m, f in zip(cols["mother"], cols["father"], strict=True) if m in row and f in row})
    x, y = [], []
    for m, f in pairs:
        xv, yv = traits[i][row[m]], traits[j][row[f]]
        if xv is not None and yv is not None:
            x.append(xv)
            y.append(yv)
    case.oracle_payload = {"x": x, "y": y, "kinds_xy": [kinds[i], kinds[j]]}
    # Every couple has its own child and no parent remates: each Mate Network is one pair.
    case.context = {"n": len(x), "G": len(x)}


def _build_po(case: Case, cols: dict[str, list[Any]]) -> None:
    """The lm data of each line and subset, rebuilt by pg's drop rules (regression.rs:338-366)."""
    stratified = bool(case.spec.get("stratified"))
    case.payload = {k: cols[k] for k in ("id", "mother", "father", "sex", "y", "stratum")} | {"stratified": stratified}
    row = {rid: i for i, rid in enumerate(cols["id"])}
    y = cols["y"]
    subsets: dict[str, list[int]] = {"pooled": list(range(len(y)))}
    if stratified:
        for s in sorted({s for s in cols["stratum"] if s is not None}):
            subsets[f"s{s}"] = [k for k, v in enumerate(cols["stratum"]) if v == s]
    fits: dict[str, dict[str, list[Any]]] = {}
    for name, members in subsets.items():
        for line in LINES:
            x_out, y_out, cluster = [], [], []
            for k in members:
                m, f = row.get(cols["mother"][k]), row.get(cols["father"][k])
                need = {"midparent": (m, f), "mother": (m,), "father": (f,)}[line]
                if any(p is None for p in need):
                    continue  # parent not in pedigree
                if y[k] is None:
                    continue  # child value missing
                parent_values = [y[p] for p in need]
                if any(v is None for v in parent_values):
                    continue  # parent value missing
                x_out.append(sum(parent_values) / len(parent_values))
                y_out.append(y[k])
                cluster.append(cols["net_via_father" if line == "father" else "net_via_mother"][k])
            key = f"{name}.{line}"
            fits[key] = {"x": x_out, "y": y_out, "cluster": cluster}
            case.context[f"{key}.G"] = len(set(cluster))
            case.context[f"{key}.n"] = len(x_out)
    case.oracle_payload = {"fits": fits}


# ---------------------------------------------------------------- pg (Python)


def pg_python(case: Case) -> dict[str, Any]:
    """The Python binding's results, under the names pg_r.R uses."""
    from pg_phenotype import Trait, assortative, correlation, pafgrs, regression

    p = case.payload
    if case.family == "tetrachoric":
        r = correlation.tetrachoric(table=p["table"])
        return {"rho": r.value, "se": r.se, "boundary": r.boundary, "n": r.n}
    if case.family == "mate":
        ped = _pedigree(p)
        traits = [Trait(_floats(t), kind=k) for t, k in zip(p["traits"], p["kinds"], strict=True)]
        res = assortative.mate_correlation(ped, traits, permutations=0)
        cell = res.cell(*p["cell"])
        r = cell.primary
        return {
            "estimator": r.estimator,
            "rho": r.value,
            "se": r.se,
            "boundary": r.boundary,
            "n": cell.n,
            "n_mate_networks": cell.n_mate_networks,
        }
    if case.family == "parent_offspring":
        stratum = _floats(p["stratum"]) if p["stratified"] else None
        res = regression.parent_offspring(_pedigree(p), Trait(_floats(p["y"]), kind="continuous"), stratum=stratum)
        out: dict[str, Any] = {}
        subsets = {"pooled": res.pooled} | {f"s{s}": sub for s, sub in (res.strata or {}).items()}
        for name, sub in subsets.items():
            for line in LINES:
                fit = getattr(sub, line)
                for f in ("slope", "se", "se_ols", "n", "n_mate_networks"):
                    out[f"{name}.{line}.{f}"] = getattr(fit, f)
        return out
    if case.family == "krebs":
        ped = {
            "id": np.array([1, 2]),
            "mother": np.array([-1, -1]),
            "father": np.array([-1, -1]),
            "twin": np.array([2, 1]),
        }
        prep = pafgrs.prepare(ped, ndegree=1)
        out = {}
        for name, s in p["scenarios"].items():
            cip = pafgrs.Cip(np.array([0.0, s["cip_age"]]), np.array([0.0, s["K"]]))
            trait = Trait(np.array([0.0, s["status"]]), kind="binary")
            r = pafgrs.score_univariate(prep, trait, age=np.array([math.nan, s["age"]]), cip=cip, h2=1.0)
            out[name] = float(r.est[0])
        return out
    if case.family == "aalen_johansen":
        entry = None if p["entry"] is None else _floats(p["entry"])
        cip = pafgrs.Cip.from_aalen_johansen(_floats(p["age"]), _floats(p["status"]), entry=entry)
        return {"ages": cip.ages.tolist(), "cip": cip.cip.tolist()}
    raise ValueError(case.family)


def _pedigree(p: dict[str, Any]) -> dict[str, np.ndarray]:
    return {k: np.array(p[k], dtype=np.int64) for k in ("id", "mother", "father", "sex")}


# ---------------------------------------------------------------- R scripts


def run_r(script: str, cases: list[dict[str, Any]]) -> dict[str, Any]:
    """Run one R script on *cases* through the JSON interface; return its output."""
    WORK.mkdir(exist_ok=True)
    stem = Path(script).stem
    inp, out = WORK / f"{stem}.in.json", WORK / f"{stem}.out.json"
    inp.write_text(json.dumps({"cases": cases}))
    out.unlink(missing_ok=True)
    proc = subprocess.run(
        ["Rscript", str(ROOT / script), str(inp), str(out)], check=False, capture_output=True, text=True
    )
    if proc.returncode != 0:
        raise RuntimeError(f"{script} exited {proc.returncode}: {proc.stderr.strip()[-2000:]}")
    return json.loads(out.read_text())


# ---------------------------------------------------------------- versions and preflight


@dataclass
class Versions:
    """What every row records: R, each oracle package, and both bindings."""

    python_pg: str
    python_core: str
    r_pg: str
    r_core: str
    r: str = ""
    oracles: dict[str, str] = field(default_factory=dict)

    def row(self, packages: list[str]) -> str:
        """The versions behind one row."""
        parts = [self.r.replace("R version ", "R ").split(" (")[0]]
        parts += [f"{p} {self.oracles[p]}" for p in packages if p in self.oracles]
        parts += [f"pg-phenotype {self.python_pg} (py)", f"pgphenotype {self.r_pg} (R)"]
        return "; ".join(parts)


def preflight() -> Versions:
    """Refuse to run unless both bindings report the pinned version (or, under the override, one core rev)."""
    import pg_phenotype

    pin = read_pin()
    out = run_r_versions()
    v = Versions(dist_version("pg-phenotype"), pg_phenotype.pg_core_rev(), out["version"], out["core_rev"], out["R"])
    v.oracles["scipy"] = scipy.__version__
    if os.environ.get(OVERRIDE_ENV):
        if v.python_core != v.r_core:
            sys.exit(f"preflight: override in use, but core revs differ: Python {v.python_core}, R {v.r_core}")
        print(f"preflight: override, both bindings at core {v.python_core} (py {v.python_pg}, R {v.r_pg})")
    else:
        bad = {k: x for k, x in (("Python", v.python_pg), ("R", v.r_pg)) if x != pin["version"]}
        if bad:
            sys.exit(
                f"preflight: pinned pg-phenotype {pin['version']}, but {bad}; run `pixi install` and `pixi run r-setup`"
            )
        if v.python_core != v.r_core:
            sys.exit(f"preflight: core revs differ: Python {v.python_core}, R {v.r_core}")
        print(f"preflight: pg-phenotype {v.python_pg} on both bindings, core {v.python_core[:12]}")
    return v


def run_r_versions() -> dict[str, str]:
    """pgphenotype_version() and pg_core_rev() of the installed R binding."""
    expr = (
        "suppressPackageStartupMessages(library(pgphenotype)); "
        "cat(pgphenotype_version(), pg_core_rev(), R.version.string, sep = '\\n')"
    )
    proc = subprocess.run(["Rscript", "-e", expr], check=False, capture_output=True, text=True)
    if proc.returncode != 0:
        sys.exit(f"preflight: the R binding did not load ({proc.stderr.strip()}); run `pixi run r-setup`")
    version, core, r = proc.stdout.strip().splitlines()
    return {"version": version, "core_rev": core, "R": r}


# ---------------------------------------------------------------- checks


def krebs_closed_form(s: dict[str, float]) -> float:
    """E[co-twin liability | status] at h2 = 1: phi(T)/K for a case, -w phi(T)/(1 - w K) for a control.

    w = CIP(age)/K with CIP linear from 0 at age 0 to K at ``cip_age``
    (pg's mixture, crates/core/src/normal.rs:234, pafgrs/cip.rs:127).
    """
    k = s["K"]
    t = norm.isf(k)
    if s["status"] == 1:
        return float(norm.pdf(t) / k)
    w = min(s["age"] / s["cip_age"], 1.0)
    return float(-w * norm.pdf(t) / (1 - w * k))


@dataclass
class Row:
    """One line of results."""

    case: str
    check: str
    reference: str = ""
    conversion: str = ""
    converted: str = ""
    pg_py: str = ""
    pg_r: str = ""
    diff: str = ""
    tol: str = ""
    status: str = ""
    id: str = ""
    source: str = ""
    note: str = ""
    versions: str = ""


def _fmt(v: Any) -> str:
    if v is None:
        return ""
    if isinstance(v, bool):
        return str(v).lower()
    if isinstance(v, float):
        return repr(v)
    if isinstance(v, list):
        return f"[{len(v)} values]"
    return str(v)


def _gap(a: Any, b: Any, relative: bool) -> float:
    """Max |a - b| (relative to |b| when *relative*) over scalars or equal-length vectors."""
    if isinstance(a, (bool, str)) or isinstance(b, (bool, str)):
        return 0.0 if a == b else math.inf
    av, bv = np.atleast_1d(np.asarray(a, dtype=float)), np.atleast_1d(np.asarray(b, dtype=float))
    if av.shape != bv.shape:
        return math.inf
    d = np.abs(av - bv)
    if relative:
        d = d / np.abs(bv)
    return float(d.max()) if d.size else 0.0


def evaluate(
    case: Case, check: dict[str, Any], py: dict[str, Any], r: dict[str, Any], oracle: dict[str, Any], versions: Versions
) -> Row:
    """One check's row."""
    expect = check.get("expect", "pass")
    row = Row(
        case=case.id,
        check=check["label"],
        id=check.get("id", ""),
        source=check.get("source", ""),
        note=check.get("note", ""),
    )
    if expect == "gap":
        row.status = "GAP"
        return row
    ref_spec = check["ref"]
    packages = []
    if "oracle" in ref_spec:
        script, key = ref_spec["oracle"].split(".", 1)
        packages = [script]
        values = oracle[script]["values"][case.id]
        if "error" in values:
            raise RuntimeError(f"{script}: {values['error']}")
        reference = values[key]
    elif "value" in ref_spec:
        reference = ref_spec["value"]
    elif "context" in ref_spec:
        reference = case.context[ref_spec["context"]]
    elif ref_spec.get("closed_form") == "krebs":
        reference = krebs_closed_form(case.payload["scenarios"][check["quantity"]])
        packages = ["scipy"]
    else:
        raise ValueError(f"unknown ref {ref_spec!r}")
    row.versions = versions.row(packages)
    row.reference = _fmt(reference)
    converted = reference
    if "conversion" in check:
        cid = check["conversion"]
        if cid not in CONVERSIONS:
            raise ValueError(f"unknown conversion {cid!r}" + (" (reserved)" if cid in RESERVED else ""))
        factor = CONVERSIONS[cid].factor(case.context)
        converted = None if reference is None else reference * factor
        row.conversion = f"{cid} x{factor:.10g}"
        row.converted = _fmt(converted)
    a, b = py[check["quantity"]], r[check["quantity"]]
    row.pg_py, row.pg_r = _fmt(a), _fmt(b)
    relative = bool(check.get("relative"))
    tol = float(check["tol"])
    row.tol = f"{tol:g}" + (" rel" if relative else "")
    notes = [row.note] if row.note else []
    bindings_ok = (a is None) == (b is None) and (a is None or _gap(a, b, False) <= BINDINGS_TOL)
    if a is None or converted is None:
        gap = math.inf
    else:
        gap = max(_gap(a, converted, relative), _gap(b, converted, relative))
    row.diff = "" if math.isinf(gap) else f"{gap:.3g}"
    within = gap <= tol
    if not bindings_ok:
        row.status = "FAIL"
        notes.append(f"bindings disagree: py {_fmt(a)}, R {_fmt(b)}")
    elif expect == "pass":
        row.status = "PASS" if within else "FAIL"
    elif expect == "diff":
        row.status = "DIFF" if not within else "FAIL"
        if within:
            notes.append("declared difference not seen: update cases.yaml and differences.md")
    elif expect == "info":
        row.status = "INFO"
    else:
        raise ValueError(f"unknown expect {expect!r}")
    row.note = "; ".join(notes)
    return row


# ---------------------------------------------------------------- report


COLUMNS = (
    ("Case", "case"),
    ("Check", "check"),
    ("Reference", "reference"),
    ("Conversion", "conversion"),
    ("Converted", "converted"),
    ("pg (Python)", "pg_py"),
    ("pg (R)", "pg_r"),
    ("Diff", "diff"),
    ("Tol", "tol"),
    ("Status", "status"),
    ("Id", "id"),
    ("Source", "source"),
    ("Note", "note"),
    ("Versions", "versions"),
)


def write_results(rows: list[Row], versions: Versions | None) -> str:
    """results.md and results.tsv; returns the totals line."""
    RESULTS.mkdir(exist_ok=True)
    counts = {s: sum(r.status == s for r in rows) for s in STATUSES}
    totals = "Totals: " + ", ".join(f"{s} {counts[s]}" for s in sorted(counts) if counts[s])
    with (RESULTS / "results.tsv").open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow([c for _, c in COLUMNS])
        w.writerows([[getattr(r, c) for _, c in COLUMNS] for r in rows])
    lines = [
        "# pg-phenotype against outside oracles and published values",
        "",
        "Generated by `run_oracles.py`. Statuses: PASS (both bindings match the reference within the tolerance), "
        "DIFF (a known difference, explained in differences.md under its id), GAP (not checked; see its id), "
        "INFO (a different estimator, shown with its measured gap and never counted as PASS or FAIL), FAIL, ERROR. "
        f"The two bindings must also agree with each other to {BINDINGS_TOL:g}. Conversions (C1, C2) are listed in "
        "conversions.py; the row shows the oracle's raw value, the factor and the converted value.",
        "",
        totals,
        "",
    ]
    if versions is not None:
        oracles = ", ".join(f"{k} {v}" for k, v in sorted(versions.oracles.items()))
        lines += [
            f"pg-phenotype {versions.python_pg} (Python, core {versions.python_core[:12]}); "
            f"pgphenotype {versions.r_pg} (R, core {versions.r_core[:12]}); {versions.r}; {oracles}.",
            "",
        ]
    shown = [c for c in COLUMNS if c[1] != "versions"]
    lines.append("| " + " | ".join(h for h, _ in shown) + " |")
    lines.append("|" + "---|" * len(shown))
    for r in rows:
        cells = [str(getattr(r, c)).replace("|", "\\|") for _, c in shown]
        lines.append("| " + " | ".join(cells) + " |")
    lines += ["", "Each row's versions (R, the oracle package, both bindings) are in `results.tsv`.", ""]
    (RESULTS / "results.md").write_text("\n".join(lines))
    return totals


# ---------------------------------------------------------------- main


ORACLE_SCRIPTS = {
    "polycor": "oracles/polycor.R",
    "psych": "oracles/psych.R",
    "sandwich": "oracles/sandwich.R",
    "survival": "oracles/survival.R",
    "cmprsk": "oracles/cmprsk.R",
}


def main() -> int:
    """Run the cases and write results/."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--install-r", action="store_true", help="install the pinned R tarball and exit")
    parser.add_argument("--only", help="comma-separated case ids")
    args = parser.parse_args()
    if args.install_r:
        install_r()
        return 0
    specs = yaml.safe_load((ROOT / "cases.yaml").read_text())["cases"]
    if args.only:
        wanted = set(args.only.split(","))
        specs = [s for s in specs if s["id"] in wanted]
    versions = preflight()
    cases = [build(s) for s in specs]
    computed = [c for c in cases if c.family != "gap"]
    rows: list[Row] = []

    py: dict[str, Any] = {}
    errors: dict[str, str] = {}
    for c in computed:
        try:
            py[c.id] = pg_python(c)
        except Exception as e:
            errors[c.id] = f"pg (Python): {type(e).__name__}: {e}"
    r_out = run_r("pg_r.R", [{"case": c.id, "family": c.family, **c.payload} for c in computed])
    for cid, values in r_out["values"].items():
        if "error" in values:
            errors.setdefault(cid, f"pg (R): {values['error']}")

    oracle: dict[str, Any] = {}
    for name, script in ORACLE_SCRIPTS.items():
        needed = [
            c
            for c in computed
            if any(ch.get("ref", {}).get("oracle", "").startswith(name + ".") for ch in c.spec["checks"])
        ]
        if not needed:
            continue
        oracle[name] = run_r(script, [{"case": c.id, "family": c.family, **c.oracle_payload} for c in needed])
        versions.r = versions.r or oracle[name]["versions"]["R"]
        versions.oracles[name] = oracle[name]["versions"][name]
    versions.r = versions.r or r_out["versions"]["R"]

    for c in cases:
        for check in c.spec["checks"]:
            try:
                if c.id in errors:
                    raise RuntimeError(errors[c.id])
                rows.append(evaluate(c, check, py.get(c.id, {}), r_out["values"].get(c.id, {}), oracle, versions))
            except Exception as e:
                rows.append(Row(case=c.id, check=check["label"], status="ERROR", note=f"{type(e).__name__}: {e}"))
    totals = write_results(rows, versions)
    print(totals)
    for r in rows:
        if r.status in ("FAIL", "ERROR"):
            print(f"{r.status}: {r.case} / {r.check}: {r.note} (diff {r.diff}, tol {r.tol})")
    return 1 if any(r.status in ("FAIL", "ERROR") for r in rows) else 0


if __name__ == "__main__":
    sys.exit(main())
