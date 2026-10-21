"""The closed list of declared conversions (ADR 0001).

A conversion turns an oracle's value into pg-phenotype's convention when the
two differ by a known factor.  Each takes the case's context (``n`` pairs,
``G`` Mate Networks) and returns the factor applied to the oracle value.  A
check naming an id not listed here is an ERROR.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Callable, Mapping


@dataclass(frozen=True)
class Conversion:
    """One declared convention gap: what it multiplies, and where pg-phenotype states it."""

    id: str
    description: str
    source: str
    factor: Callable[[Mapping[str, float]], float]


CONVERSIONS = {
    c.id: c
    for c in (
        Conversion(
            "C1",
            "SE x sqrt(n/(n-1)): pg's small-sample factor, every pair its own cluster",
            "pg-phenotype crates/core/src/correlation.rs:163",
            lambda ctx: math.sqrt(ctx["n"] / (ctx["n"] - 1)),
        ),
        Conversion(
            "C2",
            "SE x sqrt(G/(G-1)): pg's Mate Network sandwich factor",
            "pg-phenotype crates/core/src/mate_network.rs:138-142",
            lambda ctx: math.sqrt(ctx["G"] / (ctx["G"] - 1)),
        ),
    )
}
#: C3 (SD divisor n vs n-1 for polyserial) is reserved: no exact factor was
#: derived, so those rows are DIFF (differences.md D3).
RESERVED = {"C3"}
