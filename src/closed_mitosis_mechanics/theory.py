"""Exact analytical theory for confinement-driven spindle-nuclear-envelope loading.

Model. A long membrane tether (nuclear envelope, NE) obeys the Derenyi et al.
2002 free energy per unit length, plus a uniform external effective pressure
p_eff acting as a full-transmission upper bound on cortex-to-NE confinement:

    g(R)/L = pi*kappa*R^-1 + 2*pi*sigma*R + pi*p_eff*R^2

Stationarity (dg/dR=0) gives 2*p_eff*R^3 + 2*sigma*R^2 - kappa = 0. Writing
R_D = sqrt(kappa/(2*sigma)), x = R/R_D, q = p_eff*R_D/sigma reduces this to
the depressed cubic q*x^3 + x^2 - 1 = 0, and substituting back into the pull
force f(R) = pi*kappa/R + 2*pi*sigma*R + pi*p_eff*R^2 gives the closed form

    f(p_eff) = (f0(0)/4) * (3/x + x),      f0(0) = 2*pi*sqrt(2*sigma*kappa)

for whichever x solves the cubic (this also reduces to f0(0) at p_eff=0,
x=1). Setting f(p_star) = F_b and mu = F_b/f0(0) gives 3/x*+x* = 4*mu, i.e.
x*^2 - 4*mu*x* + 3 = 0, whose physical branch (x*<=1) is

    x* = 2*mu - sqrt(4*mu^2 - 3)

and the dimensionless threshold q*(mu) = (1-x*^2)/x*^3 (equivalently the
explicit polynomial-plus-radical form in q_star_exact below) gives

    p_star = q*(mu) * sqrt(2*sigma^3/kappa)

mu <= 1 semantics (2026-08-25 correction). p_star is always computed as a
well-defined, nonnegative number -- this module never raises for mu<=1.
Three physically distinct regimes are instead reported via
ThresholdResult.status / .confinement_induced_crossing:

    mu > 1  POSITIVE_PRESSURE_CROSSING  p_star>0   REQUIRED
        A positive effective pressure is required to reach F_b -- the
        ordinary confinement-induced-crossing scenario.
    mu == 1 ZERO_PRESSURE_THRESHOLD     p_star=0   NOT_REQUIRED
        F_b is reached exactly at zero additional effective pressure -- a
        knife-edge boundary, not a crossing confinement causes.
    mu < 1  BASELINE_ALREADY_SUPERTHRESHOLD p_star=0  NOT_APPLICABLE
        f0(0) already exceeds F_b with no confinement at all. There is no
        "crossing" for confinement to induce; the question this model asks
        does not apply. p_star=0 here is a semantic clamp (Tcrit_min=0),
        NOT a root of the stationarity equation -- substituting mu<1 into
        the raw q_star_exact() formula returns spurious *negative* values
        (e.g. q_star_exact(0.9) = -0.318), because that formula assumes the
        physical p_eff>=0 branch, which does not exist for mu<1.

Universal criterion. With T the (unmeasured) cortex-to-NE transmission
fraction and p_confine the actual confinement pressure, define the
effective pressure p_eff = T*p_confine and the dimensionless transmitted
load

    chi = T * p_confine * sqrt(kappa / (2*sigma^3))

Threshold compatibility holds, for any parameter values, exactly iff

    chi >= q_star(mu),   mu >= 1

(for mu<=1 the right-hand side is 0 by the clamp above, so the criterion is
trivially satisfied by any chi>=0 -- consistent with "the baseline is
already super-threshold regardless of confinement"). This is the model's
primary, parameter-free result; any single numerical p_star (e.g. the
illustrative S. pombe baseline in parameters.py) is one evaluation of it at
one F_b, not an independent finding.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

TWO_PI = 2.0 * math.pi
SQRT3_OVER_2 = math.sqrt(3.0) / 2.0

THRESHOLD_STATUS_POSITIVE_CROSSING = "POSITIVE_PRESSURE_CROSSING"
THRESHOLD_STATUS_ZERO_PRESSURE = "ZERO_PRESSURE_THRESHOLD"
THRESHOLD_STATUS_BASELINE_SUPERTHRESHOLD = "BASELINE_ALREADY_SUPERTHRESHOLD"

CROSSING_REQUIRED = "REQUIRED"
CROSSING_NOT_REQUIRED = "NOT_REQUIRED"
CROSSING_NOT_APPLICABLE = "NOT_APPLICABLE"


def f0_zero_pressure(kappa_j: float, sigma_npm: float) -> float:
    """f0(0) = 2*pi*sqrt(2*sigma*kappa), in N (kappa in J, sigma in N/m)."""
    if kappa_j <= 0.0 or sigma_npm <= 0.0:
        raise ValueError("kappa and sigma must be positive")
    return TWO_PI * math.sqrt(2.0 * sigma_npm * kappa_j)


def mu_margin(buckling_n: float, kappa_j: float, sigma_npm: float) -> float:
    """mu = F_b / f0(0), the intrinsic spindle-to-membrane force margin."""
    return buckling_n / f0_zero_pressure(kappa_j, sigma_npm)


def x_star_from_mu(mu: float) -> float:
    """Physical-branch root of x^2 - 4*mu*x + 3 = 0: x* = 2*mu - sqrt(4*mu^2-3).

    Raises ValueError for mu < sqrt(3)/2 (discriminant negative -- no real
    root at all). For sqrt(3)/2 <= mu < 1 a real root exists but is outside
    the physical p_eff>=0 branch; callers needing the mu<=1 semantic clamp
    should use p_star() below, not this function directly.
    """
    if mu < SQRT3_OVER_2:
        raise ValueError(f"mu={mu:.6f} < sqrt(3)/2: discriminant is negative")
    return 2.0 * mu - math.sqrt(4.0 * mu * mu - 3.0)


def q_star_exact(mu: float) -> float:
    """q*(mu) = (1-x*^2)/x*^3 = (4/27)*[8*mu^3 - 9*mu + (4*mu^2-3)^(3/2)], mu>=1.

    The two forms are algebraically identical (see tests); this is the raw
    formula, valid only on the physical branch mu>=1 -- it returns spurious
    *negative* values for mu<1 (e.g. q_star_exact(0.9)=-0.318). Callers that
    need a value defined for all mu should use p_star()["q_star"], which
    applies the mu<=1 semantic clamp (q*=0) instead of evaluating this
    formula out of its domain.
    """
    if mu < 1.0:
        raise ValueError(f"mu={mu:.6f} < 1: q*(mu) is only defined for mu>=1")
    return (4.0 / 27.0) * (8.0 * mu**3 - 9.0 * mu + (4.0 * mu * mu - 3.0) ** 1.5)


def dq_star_dmu(mu: float) -> float:
    """Exact dq*/dmu = 4/x*^2, mu>=1 (diverges as mu->1+)."""
    x_star = x_star_from_mu(mu)
    return 4.0 / (x_star * x_star)


def d2q_star_dmu2(mu: float) -> float:
    """Exact d^2q*/dmu^2 = 32/[x*(3-x*^2)], mu>=1 (diverges as mu->1+)."""
    x_star = x_star_from_mu(mu)
    return 32.0 / (x_star * (3.0 - x_star * x_star))


@dataclass(frozen=True)
class ThresholdResult:
    """Numerical value (always defined, always >=0) plus event semantics."""

    p_star_pa: float
    mu: float
    x_star: float
    q_star: float
    threshold_status: str
    confinement_induced_crossing: str


def p_star(kappa_j: float, sigma_npm: float, buckling_n: float) -> ThresholdResult:
    """Exact closed-form critical effective pressure, with mu<=1 handled explicitly.

    Never raises: always returns a numerical p_star_pa>=0. See module
    docstring for the mu<1 / mu==1 / mu>1 status semantics.
    """
    mu = mu_margin(buckling_n, kappa_j, sigma_npm)
    if mu < 1.0:
        return ThresholdResult(
            p_star_pa=0.0,
            mu=mu,
            x_star=1.0,
            q_star=0.0,
            threshold_status=THRESHOLD_STATUS_BASELINE_SUPERTHRESHOLD,
            confinement_induced_crossing=CROSSING_NOT_APPLICABLE,
        )
    x_star = x_star_from_mu(mu)
    q_star = (1.0 - x_star * x_star) / (x_star**3)
    r_d = math.sqrt(kappa_j / (2.0 * sigma_npm))
    p_star_pa = q_star * sigma_npm / r_d
    if mu == 1.0:
        status, crossing = THRESHOLD_STATUS_ZERO_PRESSURE, CROSSING_NOT_REQUIRED
    else:
        status, crossing = THRESHOLD_STATUS_POSITIVE_CROSSING, CROSSING_REQUIRED
    return ThresholdResult(
        p_star_pa=p_star_pa,
        mu=mu,
        x_star=x_star,
        q_star=q_star,
        threshold_status=status,
        confinement_induced_crossing=crossing,
    )


def chi_criterion(transmission_fraction: float, p_confine_pa: float, kappa_j: float, sigma_npm: float) -> float:
    """chi = T * p_confine * sqrt(kappa / (2*sigma^3)), the universal dimensionless load."""
    return transmission_fraction * p_confine_pa * math.sqrt(kappa_j / (2.0 * sigma_npm**3))


def is_threshold_compatible(chi: float, mu: float) -> bool:
    """chi >= q_star(mu), the universal parameter-free feasibility criterion.

    For mu<=1 the right-hand side is 0 by the mu<=1 clamp (see module
    docstring), so this is trivially True for any chi>=0.
    """
    if mu < 1.0:
        return True
    x_star = x_star_from_mu(mu)
    q_star = (1.0 - x_star * x_star) / (x_star**3)
    return chi >= q_star


def exact_elasticities(mu: float) -> dict[str, float]:
    """Exact (non-finite-difference) log-log elasticities of p_star at fixed x*(mu).

        E_Fb    = d ln(p*)/d ln(F_b)    = (3+x*^2)/(1-x*^2)
        E_kappa = d ln(p*)/d ln(kappa)  = -2/(1-x*^2)
        E_sigma = d ln(p*)/d ln(sigma)  = -2*x*^2/(1-x*^2)

    Exact identities (see tests): E_Fb + E_kappa + E_sigma = 1;
    E_kappa - E_sigma = -2. Only defined for mu>1 (diverges as mu->1+;
    undefined at/below mu=1, where p_star's mu<=1 clamp makes it locally
    non-differentiable in the ordinary sense).
    """
    if mu <= 1.0:
        raise ValueError(f"mu={mu:.6f} <= 1: elasticities are only defined for mu>1")
    x_star = x_star_from_mu(mu)
    denom = 1.0 - x_star * x_star
    return {
        "F_b": (3.0 + x_star * x_star) / denom,
        "kappa": -2.0 / denom,
        "sigma": -2.0 * x_star * x_star / denom,
    }


def log_log_elasticity(func, x0: float, *, step: float = 1.0e-4) -> float:
    """Central-difference d ln(func(x)) / d ln(x) at x0 (independent numeric cross-check)."""
    x_up = x0 * math.exp(step)
    x_dn = x0 * math.exp(-step)
    return (math.log(func(x_up)) - math.log(func(x_dn))) / (2.0 * step)


def static_tension_multiplier_elasticities() -> dict[str, float]:
    """M_sigma=(F_b/f0(0))^2=mu^2 is proportional to F_b^2/(kappa*sigma): exact, all mu."""
    return {"kappa": -1.0, "sigma": -1.0, "F_b": 2.0}


def critical_length(a_pn_um4: float, kappa_j: float, sigma_npm: float) -> float:
    """L* (um) solving F_b(L)=A/L^4 = f0(0): L* = (A/f0(0))^(1/4).

    A in pN*um^4, f0(0) converted to pN internally. If using this mapping,
    d ln(p*)/d ln(L) = -4 * E_Fb (chain rule on F_b=A/L^4).
    """
    f0_pn = f0_zero_pressure(kappa_j, sigma_npm) * 1.0e12
    return (a_pn_um4 / f0_pn) ** 0.25
