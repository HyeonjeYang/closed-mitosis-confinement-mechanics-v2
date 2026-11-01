# Minimal reproducible code repository for the analysis and code-generated figures.

This repository contains the analytical model and the code needed to reproduce
its numerical results, uncertainty scenarios, and code-generated figures for
external confinement during closed mitosis.

## Contents

- `src/closed_mitosis_mechanics/`: theory, parameter loading, and sensitivity routines
- `data/parameters.yaml`: model inputs and their source/status metadata
- `scripts/export_results.py`: baseline and deterministic sensitivity summaries
- `scripts/statistical_robustness.py`: independent and fitted-ratio scenario propagation
- `scripts/make_figures.py`: code-generated figures (does not include or generate the manual schematic)

Generated figures and result tables are intentionally not tracked.

## Run

```bash
python -m pip install -e .
python scripts/export_results.py
python scripts/statistical_robustness.py
python scripts/make_figures.py
```
