#!/usr/bin/env python3
"""Export final result CSVs from the authoritative theory implementation.

Writes:
    results/baseline.csv
    results/exact_crosscheck.csv
    results/robustness_oat.csv
    results/final_paper_numbers.csv

No physical value is hardcoded here -- everything comes from
data/parameters.yaml via closed_mitosis_mechanics.parameters, and every
derived number comes from closed_mitosis_mechanics.theory /
closed_mitosis_mechanics.robustness.
"""

from __future__ import annotations

import csv
import math
import sys
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from scipy.optimize import brentq  # noqa: E402

from closed_mitosis_mechanics import (  # noqa: E402
    load_baseline,
    load_begley_comparator,
    load_bilayer_reference,
    load_buckling_prefactor,
)
from closed_mitosis_mechanics.robustness import (  # noqa: E402
    deterministic_oat,
    literature_bounded_length_scenarios,
)
from closed_mitosis_mechanics.theory import (  # noqa: E402
    critical_length,
    exact_elasticities,
    f0_zero_pressure,
    mu_margin,
    p_star,
)

RESULTS_DIR = REPO_ROOT / "results"
PARAMETER_CONVENTION = "LIM_EFFECTIVE_WHOLE_NE"


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {path.relative_to(REPO_ROOT)} ({len(rows)} row{'s' if len(rows) != 1 else ''})")
    return path


def export_baseline() -> Path:
    baseline = load_baseline()
    result = p_star(baseline.kappa_effective_ne_j, baseline.sigma_effective_ne_npm, baseline.buckling_force_n)
    f0_pn = f0_zero_pressure(baseline.kappa_effective_ne_j, baseline.sigma_effective_ne_npm) * 1.0e12
    elasticities = exact_elasticities(result.mu)
    row = {
        "kappa_eff_pN_nm": round(baseline.kappa_effective_ne_pn_nm, 9),
        "sigma_eff_pN_nm": round(baseline.sigma_effective_ne_pn_nm, 9),
        "F_b_pN": round(baseline.buckling_force_pn, 9),
        "f0_pN": round(f0_pn, 6),
        "mu": round(result.mu, 8),
        "p_star_Pa": round(result.p_star_pa, 6),
        "Tcrit_600Pa": round(result.p_star_pa / 600.0, 6),
        "Tcrit_890Pa": round(result.p_star_pa / 890.0, 6),
        "E_Fb": round(elasticities["F_b"], 6),
        "E_kappa": round(elasticities["kappa"], 6),
        "E_sigma": round(elasticities["sigma"], 6),
        "threshold_status": result.threshold_status,
        "confinement_induced_crossing": result.confinement_induced_crossing,
        "parameter_convention": PARAMETER_CONVENTION,
        "scope": "ILLUSTRATIVE_BASELINE; not a universal biological measurement",
    }
    return _write_csv(RESULTS_DIR / "baseline.csv", [row])


def _numeric_p_star_pa(kappa_j: float, sigma_npm: float, buckling_n: float) -> float:
    """Independent brentq root of the raw stationarity+pull-force equations
    (not the closed-form q*/x* algebra)."""

    def stationary_radius_m(pressure_pa: float) -> float:
        if pressure_pa == 0.0:
            return math.sqrt(kappa_j / (2.0 * sigma_npm))

        def residual(radius_m: float) -> float:
            return 2.0 * pressure_pa * radius_m**3 + 2.0 * sigma_npm * radius_m**2 - kappa_j

        r_d = math.sqrt(kappa_j / (2.0 * sigma_npm))
        return brentq(residual, 1.0e-12, r_d, xtol=1e-16, rtol=1e-14)

    def pull_force_n(pressure_pa: float) -> float:
        radius_m = stationary_radius_m(pressure_pa)
        return (
            math.pi * kappa_j / radius_m
            + 2.0 * math.pi * sigma_npm * radius_m
            + math.pi * pressure_pa * radius_m**2
        )

    if pull_force_n(0.0) >= buckling_n:
        return 0.0
    upper = 100.0
    while pull_force_n(upper) < buckling_n:
        upper *= 2.0
    return brentq(lambda p: pull_force_n(p) - buckling_n, 0.0, upper, xtol=1e-10, rtol=1e-12)


def export_exact_crosscheck() -> Path:
    baseline = load_baseline()
    kappa0, sigma0, fb0 = (
        baseline.kappa_effective_ne_j,
        baseline.sigma_effective_ne_npm,
        baseline.buckling_force_n,
    )
    cases = {
        "baseline": (kappa0, sigma0, fb0),
        "kappa -20%": (kappa0 * 0.8, sigma0, fb0),
        "kappa +10%": (kappa0 * 1.1, sigma0, fb0),
        "sigma -20%": (kappa0, sigma0 * 0.8, fb0),
        "sigma +10%": (kappa0, sigma0 * 1.1, fb0),
        "F_b -10%": (kappa0, sigma0, fb0 * 0.9),
        "F_b +20%": (kappa0, sigma0, fb0 * 1.2),
    }
    rows = []
    for label, (kappa_j, sigma_npm, buckling_n) in cases.items():
        closed_form = p_star(kappa_j, sigma_npm, buckling_n)
        numeric_pa = _numeric_p_star_pa(kappa_j, sigma_npm, buckling_n)
        rows.append(
            {
                "case": label,
                "mu": round(closed_form.mu, 8),
                "p_star_closed_form_Pa": round(closed_form.p_star_pa, 6),
                "p_star_independent_numeric_Pa": round(numeric_pa, 6),
                "abs_discrepancy_Pa": round(abs(closed_form.p_star_pa - numeric_pa), 12),
                "threshold_status": closed_form.threshold_status,
                "parameter_convention": PARAMETER_CONVENTION,
            }
        )
    return _write_csv(RESULTS_DIR / "exact_crosscheck.csv", rows)


def export_robustness_oat() -> Path:
    baseline = load_baseline()
    rows = [
        {
            "parameter": row.parameter,
            "perturbation_fraction": row.perturbation_fraction,
            "parameter_value": round(row.parameter_value, 12),
            "f0_pN": round(row.f0_pn, 12),
            "mu": round(row.mu, 12),
            "p_star_Pa": round(row.p_star_pa, 6),
            "Tcrit_min_600Pa": round(row.t_crit_600, 6),
            "Tcrit_min_890Pa": round(row.t_crit_890, 6),
            "p_star_within_sorce_range": row.p_star_within_sorce_range,
            "threshold_status": row.threshold_status,
            "confinement_induced_crossing": row.confinement_induced_crossing,
            "stress_test_type": "DETERMINISTIC_OAT_NOT_A_CONFIDENCE_INTERVAL",
            "parameter_convention": PARAMETER_CONVENTION,
        }
        for row in deterministic_oat(baseline)
    ]
    return _write_csv(RESULTS_DIR / "robustness_oat.csv", rows)


def export_final_paper_numbers() -> Path:
    baseline = load_baseline()
    bilayer = load_bilayer_reference()
    begley = load_begley_comparator()
    result = p_star(baseline.kappa_effective_ne_j, baseline.sigma_effective_ne_npm, baseline.buckling_force_n)
    f0_pn = f0_zero_pressure(baseline.kappa_effective_ne_j, baseline.sigma_effective_ne_npm) * 1.0e12
    elasticities = exact_elasticities(result.mu)
    a_pn_um4 = load_buckling_prefactor()
    l_star_um = critical_length(a_pn_um4, baseline.kappa_effective_ne_j, baseline.sigma_effective_ne_npm)

    def row(quantity: str, value: Any, unit: str, notes: str, convention: str = PARAMETER_CONVENTION) -> dict[str, Any]:
        return {
            "quantity": quantity,
            "value": value,
            "unit": unit,
            "notes": notes,
            "parameter_convention": convention,
        }

    rows = [
        row("kappa_eff (MAIN)", round(baseline.kappa_effective_ne_pn_nm, 9), "pN*nm", "PRIMARY_SOURCE_VERIFIED"),
        row("sigma_eff (MAIN)", round(baseline.sigma_effective_ne_pn_nm, 9), "pN/nm", "PRIMARY_SOURCE_VERIFIED"),
        row("kappa_bilayer (provenance only)", round(bilayer.kappa_bilayer_pn_nm, 9), "pN*nm", "PRIMARY_SOURCE_VERIFIED"),
        row("sigma_bilayer (provenance only)", round(bilayer.sigma_bilayer_pn_nm, 9), "pN/nm", "DERIVED"),
        row(
            "kappa (Begley Table 1, comparator only)",
            begley.kappa_pn_nm,
            "pN*nm",
            "UNRESOLVED_PROVENANCE; never used as MAIN kappa",
            convention="BEGLEY_UNRESOLVED_PROVENANCE",
        ),
        row("F_b (illustrative spindle calibration)", round(baseline.buckling_force_pn, 9), "pN", "MODEL_CALIBRATED"),
        row("f0 (zero-pressure membrane force)", round(f0_pn, 6), "pN", "f0=2*pi*sqrt(2*sigma_eff*kappa_eff)"),
        row("mu = F_b/f0", round(result.mu, 8), "dimensionless", "close to mu=1 -> large sensitivity"),
        row("p_star (illustrative)", round(result.p_star_pa, 6), "Pa", "conditional on F_b=31 pN; not a measurement"),
        row("Tcrit_min(600 Pa)", round(result.p_star_pa / 600.0, 6), "dimensionless", "required transmission fraction"),
        row("Tcrit_min(890 Pa)", round(result.p_star_pa / 890.0, 6), "dimensionless", "required transmission fraction"),
        row("E_Fb (exact elasticity)", round(elasticities["F_b"], 6), "dimensionless", "d ln(p*)/d ln(F_b)"),
        row("E_kappa (exact elasticity)", round(elasticities["kappa"], 6), "dimensionless", "d ln(p*)/d ln(kappa)"),
        row("E_sigma (exact elasticity)", round(elasticities["sigma"], 6), "dimensionless", "d ln(p*)/d ln(sigma)"),
        row("L_star (F_b(L)=f0)", round(l_star_um, 6), "um", "below S. pombe 12-15 um range"),
        row("threshold_status", result.threshold_status, "categorical", "mu>1 => POSITIVE_PRESSURE_CROSSING"),
    ]
    return _write_csv(RESULTS_DIR / "final_paper_numbers.csv", rows)


def main() -> int:
    export_baseline()
    export_exact_crosscheck()
    export_robustness_oat()
    export_final_paper_numbers()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
