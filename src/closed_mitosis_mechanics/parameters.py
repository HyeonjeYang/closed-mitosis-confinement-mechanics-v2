"""Parameter loading and unit conversion.

Reads only from data/parameters.yaml -- no physical value is hardcoded
anywhere else in this package. Names are explicit about which convention a
quantity belongs to (kappa_bilayer vs kappa_effective_NE) so the two
per-bilayer/whole-NE conventions cannot be silently mixed.

Guard against accidentally selecting the UNRESOLVED_PROVENANCE Begley
kappa=40 value as the main bending rigidity: `load_baseline()` reads its
kappa/sigma from exactly one fixed pair of YAML paths
(nuclear_envelope.kappa_effective_NE / sigma_effective_NE), which are never
the Begley entry. Reading the Begley value at all requires calling the
separate, explicitly named `load_begley_comparator()`.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

from .provenance import BEGLEY_2026_DOI, LIM_2007_DOI, LIM_HUBER_2009_DOI, WARD_2014_DOI, Status

PACKAGE_DIR = Path(__file__).resolve().parent
REPO_ROOT = PACKAGE_DIR.parent.parent
DEFAULT_PARAMS_PATH = REPO_ROOT / "data" / "parameters.yaml"

# SI base-unit conversion factors. Only unit *conversion* constants live
# here (universal, not literature values).
_UNIT_TO_SI: dict[str, tuple[str, float]] = {
    "Pa": ("Pa", 1.0),
    "mN/m": ("N/m", 1.0e-3),
    "pN/nm": ("N/m", 1.0e-3),  # 1 pN/nm = 1e-12 N / 1e-9 m = 1e-3 N/m
    "nm": ("m", 1.0e-9),
    "um": ("m", 1.0e-6),
    "pN": ("N", 1.0e-12),
    "pN*nm": ("J", 1.0e-12 * 1.0e-9),  # 1 pN*nm = 1e-21 J
    "pN*um^4": ("pN*um^4", 1.0),  # left in mixed units; only used via theory.critical_length
    "um^-2": ("m^-2", 1.0e12),
    "dimensionless": ("dimensionless", 1.0),
}


class ParameterError(ValueError):
    """A required parameter is missing, has an unexpected unit, or fails a
    status/value sanity check."""


def to_si(value: float, unit: str) -> tuple[float, str]:
    """Convert a scalar value to SI (or the documented mixed unit for pN*um^4)."""
    if unit not in _UNIT_TO_SI:
        raise ParameterError(f"no SI conversion registered for unit {unit!r}")
    si_unit, factor = _UNIT_TO_SI[unit]
    return value * factor, si_unit


def load_params(path: Path | None = None) -> dict[str, Any]:
    params_path = path or DEFAULT_PARAMS_PATH
    with open(params_path, "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def get_entry(params: dict[str, Any], dotted_path: str) -> Any:
    node: Any = params
    for key in dotted_path.split("."):
        if not isinstance(node, dict) or key not in node:
            return None
        node = node[key]
    return node


def _required_scalar(
    params: dict[str, Any], dotted_path: str, expected_unit: str, allowed_status: tuple[str, ...]
) -> tuple[float, str]:
    entry = get_entry(params, dotted_path)
    if not isinstance(entry, dict) or entry.get("value") is None:
        raise ParameterError(f"required parameter {dotted_path!r} is missing")
    if entry.get("unit") != expected_unit:
        raise ParameterError(f"{dotted_path!r} must use unit {expected_unit!r}, got {entry.get('unit')!r}")
    if entry.get("status") not in allowed_status:
        raise ParameterError(
            f"{dotted_path!r} has status {entry.get('status')!r}, not in allowed set {allowed_status}"
        )
    if not entry.get("source_doi"):
        raise ParameterError(f"{dotted_path!r} has no source_doi")
    value_si, si_unit = to_si(float(entry["value"]), expected_unit)
    return value_si, si_unit


# Fixed YAML paths for the MAIN (Lim effective whole-NE) baseline. These are
# the only paths load_baseline() ever reads kappa/sigma from.
KAPPA_EFFECTIVE_NE_PATH = "nuclear_envelope.kappa_effective_NE"
SIGMA_EFFECTIVE_NE_PATH = "nuclear_envelope.sigma_effective_NE"
KAPPA_BILAYER_PATH = "nuclear_envelope.kappa_bilayer"
SIGMA_BILAYER_PATH = "nuclear_envelope.sigma_bilayer"
BEGLEY_KAPPA_PATH = "nuclear_envelope.begley_kappa_unresolved"
BEGLEY_SIGMA_PATH = "nuclear_envelope.begley_sigma"
BUCKLING_FORCE_PATH = "spindle.buckling_force_illustrative"
BUCKLING_PREFACTOR_PATH = "spindle.buckling_prefactor_A"
CONFINEMENT_RANGE_PATH = "confinement.epithelial_rounding_pressure"


@dataclass(frozen=True)
class Baseline:
    """The MAIN illustrative baseline: Lim effective whole-NE kappa/sigma
    paired with the adopted spindle buckling-force threshold. Every field
    here is ILLUSTRATIVE_BASELINE, not a universal biological measurement
    (see the source and status fields in data/parameters.yaml)."""

    kappa_effective_ne_j: float
    sigma_effective_ne_npm: float
    buckling_force_n: float
    kappa_effective_ne_pn_nm: float
    sigma_effective_ne_pn_nm: float
    buckling_force_pn: float


def load_baseline(params: dict[str, Any] | None = None) -> Baseline:
    """Load the MAIN baseline (kappa_effective_NE=400 pN*nm, sigma_effective_NE
    =0.026 pN/nm, F_b=31 pN). Never resolves the Begley kappa=40 value --
    that requires the separate load_begley_comparator()."""
    params = params if params is not None else load_params()
    kappa_j, _ = _required_scalar(
        params, KAPPA_EFFECTIVE_NE_PATH, "pN*nm", (Status.PRIMARY_SOURCE_VERIFIED,)
    )
    sigma_npm, _ = _required_scalar(
        params, SIGMA_EFFECTIVE_NE_PATH, "pN/nm", (Status.PRIMARY_SOURCE_VERIFIED,)
    )
    buckling_n, _ = _required_scalar(
        params, BUCKLING_FORCE_PATH, "pN", (Status.MODEL_CALIBRATED,)
    )
    return Baseline(
        kappa_effective_ne_j=kappa_j,
        sigma_effective_ne_npm=sigma_npm,
        buckling_force_n=buckling_n,
        kappa_effective_ne_pn_nm=kappa_j / 1.0e-21,
        sigma_effective_ne_pn_nm=sigma_npm / 1.0e-3,
        buckling_force_pn=buckling_n / 1.0e-12,
    )


@dataclass(frozen=True)
class BilayerReference:
    """Per-bilayer provenance quantities: PROVENANCE_ONLY, never substituted
    for the effective whole-NE pair in a main calculation."""

    kappa_bilayer_j: float
    sigma_bilayer_npm: float
    kappa_bilayer_pn_nm: float
    sigma_bilayer_pn_nm: float


def load_bilayer_reference(params: dict[str, Any] | None = None) -> BilayerReference:
    params = params if params is not None else load_params()
    kappa_j, _ = _required_scalar(params, KAPPA_BILAYER_PATH, "pN*nm", (Status.PRIMARY_SOURCE_VERIFIED,))
    sigma_npm, _ = _required_scalar(params, SIGMA_BILAYER_PATH, "pN/nm", (Status.DERIVED,))
    return BilayerReference(
        kappa_bilayer_j=kappa_j,
        sigma_bilayer_npm=sigma_npm,
        kappa_bilayer_pn_nm=kappa_j / 1.0e-21,
        sigma_bilayer_pn_nm=sigma_npm / 1.0e-3,
    )


@dataclass(frozen=True)
class BegleyComparator:
    """Begley et al. 2026 Table 1's own (kappa, sigma) pair: UNRESOLVED_PROVENANCE.
    Reproduces Begley's own reported f0=6.4 pN exactly when paired together,
    but is never selectable as the MAIN kappa/sigma by load_baseline()."""

    kappa_pn_nm: float
    sigma_pn_nm: float
    kappa_j: float
    sigma_npm: float
    status: str


def load_begley_comparator(params: dict[str, Any] | None = None) -> BegleyComparator:
    params = params if params is not None else load_params()
    entry = get_entry(params, BEGLEY_KAPPA_PATH)
    if not isinstance(entry, dict) or entry.get("status") != Status.UNRESOLVED_PROVENANCE:
        raise ParameterError(
            f"{BEGLEY_KAPPA_PATH!r} must carry status={Status.UNRESOLVED_PROVENANCE!r}"
        )
    kappa_j, _ = to_si(float(entry["value"]), "pN*nm")
    sigma_entry = get_entry(params, BEGLEY_SIGMA_PATH)
    sigma_npm, _ = to_si(float(sigma_entry["value"]), "pN/nm")
    return BegleyComparator(
        kappa_pn_nm=float(entry["value"]),
        sigma_pn_nm=float(sigma_entry["value"]),
        kappa_j=kappa_j,
        sigma_npm=sigma_npm,
        status=entry["status"],
    )


def load_buckling_prefactor(params: dict[str, Any] | None = None) -> float:
    """A in F_b(L)=A/L^4, pN*um^4 (Begley 2026's extraction from Ward 2014 Fig. 4E)."""
    params = params if params is not None else load_params()
    entry = get_entry(params, BUCKLING_PREFACTOR_PATH)
    if not isinstance(entry, dict) or entry.get("value") is None:
        raise ParameterError(f"required parameter {BUCKLING_PREFACTOR_PATH!r} is missing")
    return float(entry["value"])


def load_confinement_range_pa(params: dict[str, Any] | None = None) -> tuple[float, float]:
    """CONTEXT_ONLY Sorce et al. 2015 epithelial-rounding pressure range, Pa."""
    params = params if params is not None else load_params()
    entry = get_entry(params, CONFINEMENT_RANGE_PATH)
    if not isinstance(entry, dict) or entry.get("value") is None:
        raise ParameterError(f"required parameter {CONFINEMENT_RANGE_PATH!r} is missing")
    low, high = entry["value"]
    return float(low), float(high)


__all__ = [
    "Baseline",
    "BegleyComparator",
    "BilayerReference",
    "ParameterError",
    "get_entry",
    "load_baseline",
    "load_begley_comparator",
    "load_bilayer_reference",
    "load_buckling_prefactor",
    "load_confinement_range_pa",
    "load_params",
    "to_si",
    "BEGLEY_2026_DOI",
    "LIM_2007_DOI",
    "LIM_HUBER_2009_DOI",
    "WARD_2014_DOI",
]
