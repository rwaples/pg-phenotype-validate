"""Write the seeded inputs in inputs/*.tsv.gz.

Each generator is seeded, and gzip writes no timestamp, so a rerun is byte for
byte the same.  ``--verify`` regenerates into a temporary directory and
compares each file's sha256 with the one cases.yaml pins.

    python make_inputs.py            # write inputs/ and print each sha256
    python make_inputs.py --verify   # exit 1 if any file differs from cases.yaml
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import math
import sys
import tempfile
from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parent


def _cell(v: object) -> str:
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return ""
    if isinstance(v, float):
        return repr(v)
    return str(v)


def _write(path: Path, columns: dict[str, list]) -> None:
    names = list(columns)
    n = len(columns[names[0]])
    lines = ["\t".join(names)]
    lines += ["\t".join(_cell(columns[c][i]) for c in names) for i in range(n)]
    data = ("\n".join(lines) + "\n").encode()
    with path.open("wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", mtime=0, filename="") as f:
        f.write(data)


def _cut(liability: np.ndarray, thresholds: list[float]) -> list[float]:
    return [float(v) for v in np.searchsorted(thresholds, liability)]


def _mate_pairs(n: int, rho: float, seed: int, t0, t1) -> dict[str, list]:
    """N one-child couples: every Mate Network is one pair, so G = n.

    Each parent has a liability L; mothers' and fathers' liabilities have
    correlation *rho*.  ``t0(L, rng)`` and ``t1(L, rng)`` give a parent's two
    trait values.  Rows: mothers, then fathers, then children (traits
    missing).
    """
    rng = np.random.default_rng(seed)
    lm = rng.standard_normal(n)
    lf = rho * lm + math.sqrt(1 - rho**2) * rng.standard_normal(n)
    parents = np.concatenate([lm, lf])
    v0, v1 = t0(parents, rng), t1(parents, rng)
    return {
        "id": list(range(1, 3 * n + 1)),
        "mother": [-1] * (2 * n) + list(range(1, n + 1)),
        "father": [-1] * (2 * n) + list(range(n + 1, 2 * n + 1)),
        "sex": [0] * n + [1] * n + [int(s) for s in rng.integers(0, 2, n)],
        "t0": list(v0) + [None] * n,
        "t1": list(v1) + [None] * n,
    }


def _noisy(liability: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    return 0.8 * liability + 0.6 * rng.standard_normal(liability.size)


def mate_ordinal_3x4() -> dict[str, list]:
    """Mothers' 3-level t0 against fathers' 4-level t1 (cell 0, 1)."""
    return _mate_pairs(
        600,
        0.35,
        101,
        lambda L, rng: _cut(L, [-0.5, 0.7]),
        lambda L, rng: _cut(_noisy(L, rng), [-0.8, 0.0, 0.9]),
    )


def mate_ordinal_4x4() -> dict[str, list]:
    """One 4-level trait, mothers against fathers (cell 0, 0)."""
    return _mate_pairs(
        800,
        0.25,
        102,
        lambda L, rng: _cut(L, [-1.0, 0.1, 1.2]),
        lambda L, rng: _cut(_noisy(L, rng), [0.0]),
    )


def mate_polyserial() -> dict[str, list]:
    """Mothers' continuous t0 against fathers' 4-level t1 (cell 0, 1)."""
    return _mate_pairs(
        500,
        0.4,
        103,
        lambda L, rng: [float(v) for v in 10 + 2 * L],
        lambda L, rng: _cut(_noisy(L, rng), [-0.6, 0.3, 1.1]),
    )


def mate_biserial() -> dict[str, list]:
    """Mothers' continuous t0 against fathers' binary t1 (cell 0, 1)."""
    return _mate_pairs(
        500,
        0.3,
        104,
        lambda L, rng: [float(v) for v in L + 0.3 * rng.standard_normal(L.size)],
        lambda L, rng: _cut(_noisy(L, rng), [0.8]),
    )


def po_remating() -> dict[str, list]:
    """Couples with 1-4 children, 40 remating fathers, single parents and missing values.

    Each child's cluster label comes from the construction, not from
    pg-phenotype: mothers mate once, so a couple's Mate Network is its
    father's, and a parent in no Mating Pair is a network of its own.
    ``net_via_mother`` is the cluster of the midparent and mother lines,
    ``net_via_father`` of the father line (empty where that parent is not a
    row).
    """
    rng = np.random.default_rng(201)
    n_couples, n_remate, n_single = 300, 40, 15
    n_fathers = n_couples - n_remate
    names = ("id", "mother", "father", "sex", "y", "stratum", "net_via_mother", "net_via_father")
    rows: dict[str, list] = {k: [] for k in names}
    net_of: dict[int, str] = {}
    effect: dict[str, float] = {}

    def add(mother: int, father: int, sex: int, y: float) -> int:
        rid = len(rows["id"]) + 1
        for key, value in zip(
            names,
            (rid, mother, father, sex, y, int(rng.integers(0, 2)), net_of.get(mother), net_of.get(father)),
            strict=True,
        ):
            rows[key].append(value)
        return rid

    def parent(sex: int, net: str | None = None) -> int:
        rid = add(-1, -1, sex, float(rng.standard_normal()))
        net_of[rid] = net or f"S{rid}"
        effect.setdefault(net_of[rid], 0.7 * float(rng.standard_normal()))
        return rid

    fathers = [parent(1, f"N{len(rows['id']) + 1}") for _ in range(n_fathers)]
    remating = rng.choice(n_fathers, n_remate, replace=False)
    couples = [(parent(0, net_of[f]), f) for f in fathers + [fathers[i] for i in remating]]
    single_mothers = [parent(0) for _ in range(n_single)]
    single_fathers = [parent(1) for _ in range(n_single)]
    for p in rng.choice([m for m, _ in couples] + fathers, 12, replace=False):
        rows["y"][p - 1] = math.nan  # a parent value missing

    def child(mother: int, father: int) -> None:
        net = net_of[mother] if mother > 0 else net_of[father]
        observed = [rows["y"][p - 1] for p in (mother, father) if p > 0]
        mean = 0.35 * sum(v for v in observed if not math.isnan(v))
        y = mean + effect[net] + float(rng.standard_normal())
        add(mother, father, int(rng.integers(0, 2)), math.nan if rng.random() < 0.05 else y)

    for m, f in couples:
        for _ in range(int(rng.integers(1, 5))):
            child(m, f)
    for m in single_mothers:
        for _ in range(int(rng.integers(1, 4))):
            child(m, -1)
    for f in single_fathers:
        for _ in range(int(rng.integers(1, 4))):
            child(-1, f)
    return rows


def _aj(n: int, seed: int, *, deaths: bool, entry: bool) -> dict[str, list]:
    """Integer ages, so cases, deaths and censorings tie at one age."""
    rng = np.random.default_rng(seed)
    onset = rng.exponential(120, n)
    death = rng.exponential(90, n) if deaths else np.full(n, np.inf)
    censor = rng.uniform(5, 90, n)
    t = np.minimum.reduce([onset, death, censor])
    age = np.ceil(t)
    status = np.where(t == onset, 1, np.where(t == death, 2, 0))
    out: dict[str, list] = {}
    if entry:
        start = np.floor(rng.uniform(0, 1, n) * np.minimum(age - 1, 40))
        start[rng.random(n) < 0.4] = 0.0
        idle = np.flatnonzero(status == 0)[:5]  # censored at entry: at risk nowhere
        start[idle] = age[idle]
        assert not ((start >= age) & (status != 0)).any(), "a case or death needs entry < age"
        out["entry"] = [float(v) for v in start]
    out["age"] = [float(v) for v in age]
    out["status"] = [int(v) for v in status]
    return out


def aj_ties() -> dict[str, list]:
    """Cases, deaths and censorings tied at integer ages; no late entry."""
    return _aj(400, 301, deaths=True, entry=False)


def aj_no_deaths() -> dict[str, list]:
    """Cases and censorings only; no late entry."""
    return _aj(300, 302, deaths=False, entry=False)


def aj_late_entry() -> dict[str, list]:
    """Late entry, with ties and some people censored at entry."""
    return _aj(400, 303, deaths=True, entry=True)


GENERATORS = {
    f.__name__: f
    for f in (
        mate_ordinal_3x4,
        mate_ordinal_4x4,
        mate_polyserial,
        mate_biserial,
        po_remating,
        aj_ties,
        aj_no_deaths,
        aj_late_entry,
    )
}


def sha256(path: Path) -> str:
    """The file's sha256, hex."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_all(out: Path) -> dict[str, str]:
    """Write every input into *out*; return name -> sha256."""
    out.mkdir(parents=True, exist_ok=True)
    hashes = {}
    for name, generate in GENERATORS.items():
        path = out / f"{name}.tsv.gz"
        _write(path, generate())
        hashes[name] = sha256(path)
    return hashes


def pinned() -> dict[str, str]:
    """Input name -> sha256, as cases.yaml pins them."""
    cases = yaml.safe_load((ROOT / "cases.yaml").read_text())["cases"]
    return {Path(c["input"]["file"]).name.removesuffix(".tsv.gz"): c["input"]["sha256"] for c in cases if "input" in c}


def main() -> int:
    """Write or verify the inputs."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--verify", action="store_true", help="compare a fresh run with cases.yaml")
    args = parser.parse_args()
    if not args.verify:
        for name, digest in write_all(ROOT / "inputs").items():
            print(f"{digest}  inputs/{name}.tsv.gz")
        return 0
    with tempfile.TemporaryDirectory() as tmp:
        fresh = write_all(Path(tmp))
    want = pinned()
    bad = 0
    for name in sorted(set(fresh) | set(want)):
        ok = fresh.get(name) == want.get(name)
        bad += not ok
        print(f"{'ok  ' if ok else 'DIFF'} {name}: fresh {fresh.get(name)} pinned {want.get(name)}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
