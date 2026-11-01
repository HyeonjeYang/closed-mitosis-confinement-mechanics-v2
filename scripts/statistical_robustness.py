"""Compare independent and ratio-fixed uncertainty scenarios.

This is not an inferential analysis: the repository contains no sampling
distribution for kappa or F_b. Independent and fitted-ratio-preserving
uniform perturbations are used only as explicitly labeled scenarios.
"""

from __future__ import annotations

import csv
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from closed_mitosis_mechanics.parameters import load_baseline

OUTPUT = ROOT / "results" / "robustness_uniform_scenarios.csv"
SEED = 20261101
N_DRAWS = 1_000_000
BOUNDS = (0.05, 0.10, 0.20)
PRESSURE_CONTEXT_PA = 890.0


def summarize(bound: float, rng: np.random.Generator, baseline) -> list[dict[str, float]]:
    kappa = baseline.kappa_effective_ne_j
    sigma = baseline.sigma_effective_ne_npm
    fb = baseline.buckling_force_n
    k = kappa * rng.uniform(1 - bound, 1 + bound, N_DRAWS)
    force = fb * rng.uniform(1 - bound, 1 + bound, N_DRAWS)
    s_independent = sigma * rng.uniform(1 - bound, 1 + bound, N_DRAWS)
    s_ratio_fixed = (sigma / kappa) * k

    rows = []
    for scenario, sigma_draws in (
        ("fixed_fitted_sigma_over_kappa_fb", s_ratio_fixed),
        ("independent_kappa_sigma_fb", s_independent),
    ):
        f0 = 2 * np.pi * np.sqrt(2 * sigma_draws * k)
        mu = force / f0
        positive = mu > 1
        x = np.ones(N_DRAWS)
        x[positive] = 2 * mu[positive] - np.sqrt(4 * mu[positive] ** 2 - 3)
        q = np.zeros(N_DRAWS)
        q[positive] = (1 - x[positive] ** 2) / x[positive] ** 3
        p_star = q * np.sqrt(2 * sigma_draws**3 / k)

        p = float(positive.mean())
        mc_se = float(np.sqrt(p * (1 - p) / N_DRAWS))
        positive_quantiles = np.quantile(p_star[positive], [0.025, 0.5, 0.975])
        all_quantiles = np.quantile(p_star, [0.025, 0.5, 0.975])
        rows.append({
            "scenario": scenario,
            "uniform_half_width_fraction": bound,
            "draws": N_DRAWS,
            "positive_crossing_fraction": p,
            "crossing_fraction_mc95_low": max(0.0, p - 1.96 * mc_se),
            "crossing_fraction_mc95_high": min(1.0, p + 1.96 * mc_se),
            "p_star_all_q025_Pa": float(all_quantiles[0]),
            "p_star_all_median_Pa": float(all_quantiles[1]),
            "p_star_all_q975_Pa": float(all_quantiles[2]),
            "p_star_positive_q025_Pa": float(positive_quantiles[0]),
            "p_star_positive_median_Pa": float(positive_quantiles[1]),
            "p_star_positive_q975_Pa": float(positive_quantiles[2]),
            "fraction_p_star_le_890Pa": float(np.mean(p_star <= PRESSURE_CONTEXT_PA)),
            "seed": SEED,
            "distribution_assumption": scenario,
        })
    return rows


def main() -> None:
    baseline = load_baseline()
    rng = np.random.default_rng(SEED)
    rows = [row for bound in BOUNDS for row in summarize(bound, rng, baseline)]
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {OUTPUT.relative_to(ROOT)} ({N_DRAWS:,} draws per scenario, seed {SEED}).")


if __name__ == "__main__":
    main()
