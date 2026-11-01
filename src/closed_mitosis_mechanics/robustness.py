"""Deterministic one-at-a-time (OAT) robustness stress tests.

These are stress tests, not confidence intervals: no probability
distribution or joint/correlated perturbation is constructed anywhere in
this module. Single source of truth for the OAT grid used by the analysis and figure scripts.
"""

from __future__ import annotations

from dataclasses import dataclass

from .parameters import Baseline
from .theory import f0_zero_pressure, mu_margin, p_star

PERTURBATIONS: tuple[float, ...] = (-0.20, -0.10, -0.05, 0.0, 0.05, 0.10, 0.20)
REFERENCE_QUERY_PRESSURES_PA: tuple[float, float] = (600.0, 890.0)
SORCE_RANGE_PA: tuple[float, float] = (0.0, 890.0)


@dataclass(frozen=True)
class OatRow:
    parameter: str
    perturbation_fraction: float
    parameter_value: float
    f0_pn: float
    mu: float
    p_star_pa: float
    t_crit_600: float
    t_crit_890: float
    p_star_within_sorce_range: bool
    threshold_status: str
    confinement_induced_crossing: str


def _oat_row(parameter: str, fraction: float, kappa_j: float, sigma_npm: float, buckling_n: float) -> OatRow:
    result = p_star(kappa_j, sigma_npm, buckling_n)
    f0_pn = f0_zero_pressure(kappa_j, sigma_npm) * 1.0e12
    varied_value = {
        "kappa": kappa_j / 1.0e-21,
        "sigma": sigma_npm * 1.0e3,
        "F_b": buckling_n * 1.0e12,
    }[parameter]
    return OatRow(
        parameter=parameter,
        perturbation_fraction=fraction,
        parameter_value=varied_value,
        f0_pn=f0_pn,
        mu=mu_margin(buckling_n, kappa_j, sigma_npm),
        p_star_pa=result.p_star_pa,
        t_crit_600=result.p_star_pa / REFERENCE_QUERY_PRESSURES_PA[0],
        t_crit_890=result.p_star_pa / REFERENCE_QUERY_PRESSURES_PA[1],
        p_star_within_sorce_range=SORCE_RANGE_PA[0] <= result.p_star_pa <= SORCE_RANGE_PA[1],
        threshold_status=result.threshold_status,
        confinement_induced_crossing=result.confinement_induced_crossing,
    )


def deterministic_oat(baseline: Baseline) -> list[OatRow]:
    """Perturb kappa, sigma, F_b independently by +/-5/10/20%, one at a time,
    holding the other two at their ILLUSTRATIVE_BASELINE reference value."""
    kappa0, sigma0, fb0 = (
        baseline.kappa_effective_ne_j,
        baseline.sigma_effective_ne_npm,
        baseline.buckling_force_n,
    )
    rows: list[OatRow] = []
    for parameter in ("kappa", "sigma", "F_b"):
        for fraction in PERTURBATIONS:
            kappa_j = kappa0 * (1.0 + fraction) if parameter == "kappa" else kappa0
            sigma_npm = sigma0 * (1.0 + fraction) if parameter == "sigma" else sigma0
            buckling_n = fb0 * (1.0 + fraction) if parameter == "F_b" else fb0
            rows.append(_oat_row(parameter, fraction, kappa_j, sigma_npm, buckling_n))
    return rows


@dataclass(frozen=True)
class LengthScenarioRow:
    scenario: str
    length_um: float
    f_b_pn: float
    f0_pn: float
    mu: float
    p_star_pa: float
    t_crit_890: float
    p_star_within_sorce_range: bool
    threshold_status: str


def literature_bounded_length_scenarios(baseline: Baseline, a_pn_um4: float) -> list[LengthScenarioRow]:
    """F_b(L)=A/L^4 evaluated at Begley et al. 2026's own stated buckling-onset
    length window (4-5 um), a literature-motivated scenario pair, not a
    confidence interval."""
    kappa_j, sigma_npm = baseline.kappa_effective_ne_j, baseline.sigma_effective_ne_npm
    f0_pn = f0_zero_pressure(kappa_j, sigma_npm) * 1.0e12
    rows: list[LengthScenarioRow] = []
    for label, length_um in (("ONSET_WINDOW_HIGH_END", 5.0), ("ONSET_WINDOW_LOW_END", 4.0)):
        f_b_pn = a_pn_um4 / (length_um**4)
        result = p_star(kappa_j, sigma_npm, f_b_pn * 1.0e-12)
        rows.append(
            LengthScenarioRow(
                scenario=label,
                length_um=length_um,
                f_b_pn=f_b_pn,
                f0_pn=f0_pn,
                mu=result.mu,
                p_star_pa=result.p_star_pa,
                t_crit_890=result.p_star_pa / REFERENCE_QUERY_PRESSURES_PA[1],
                p_star_within_sorce_range=SORCE_RANGE_PA[0] <= result.p_star_pa <= SORCE_RANGE_PA[1],
                threshold_status=result.threshold_status,
            )
        )
    return rows


__all__ = [
    "OatRow",
    "LengthScenarioRow",
    "PERTURBATIONS",
    "REFERENCE_QUERY_PRESSURES_PA",
    "SORCE_RANGE_PA",
    "deterministic_oat",
    "literature_bounded_length_scenarios",
]
