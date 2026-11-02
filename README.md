# Temporal fidelity in wearable sweat sensing: code and data

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22932238.svg)](https://doi.org/10.5281/zenodo.22932238)

Code and machine-readable outputs for the Perspective *From Sensitivity to
Temporal Fidelity in Wearable Sweat Sensing* (submitted to ACS Sensors).

## Scope

These materials support a Perspective, that is, a peer-reviewed opinion and
review article. They contain no new experimental data. Everything here is
either a model output computed from the equations stated in the article, or a
set of coordinates digitized from figures in already published papers. Model
parameters are declared for illustration and are not measurements of any
device, and none of the numbers here should be read as an experimental
characterization of a sensor.

## Reproducing the analyses

```bash
pip install -r requirements.txt
cd code/python
python validate_transport_fidelity_models.py   # numerical checks
python run_evidence_closure_analysis.py        # canonical grids (1,080 and 1,620 runs)
python compare_figure2_strategies.py           # dense boundary grid behind Fig. 2c
python run_waveform_topology_analysis.py       # empirical stress test (1,295 runs)
python run_round2_analyses.py                  # envelopes at two sensor constants,
                                               # volume budget, time-varying flux,
                                               # fill-measure-clear and window results
python run_analytic_laplace.py                 # Laplace-domain solutions: exact periodic
                                               # envelopes, ramp and pulse rules, convergence
```

Scripts resolve paths relative to the repository root, so they run unchanged
after cloning. Figure scripts (`plot_fig1*.py`, `build_figure2_placement_ready.py`,
`rebuild_fig2b_contrast.py`, `build_fig2b_panel_variants.py`, `plot_fig3_topology.py`,
`plot_fig4_*.py`, `export_fig4_curves.py`) regenerate the plotted curves; the final
figures were assembled from these outputs.

## Layout

```
code/python/          models, analyses, figure scripts
outputs/
  revision_2026-08-31_figure2_data_package/    released grids and empirical results
  figure_working/figure2_strategy_comparison_2026-09-09/data/   dense boundary grid
  revision_2026-09-22/round2_analyses/         second-round analyses
  revision_2026-10-27/analytic/                Laplace-domain solutions (minor revision)
  figures_final/figure2_placement_ready/data/  per-panel curve values
  revision_2026-09-1*/ and 2026-09-2*/         per-figure values and parameters
```

## Key files

| File | Contents |
|---|---|
| `transport_fidelity_models.py` | plug-flow, well-mixed and gamma RTD transport, first-order sensor, fill--measure--clear and chronological-window models |
| `01_canonical_primary_grid_1080.csv` | canonical grid, 3 waveforms x 5 timescales x 3 topologies x 6 transport times x 4 sensor constants |
| `02_canonical_extended_grid_1620.csv` | same grid extended to transport times of 0.05, 0.10 and 0.25 min |
| `run_analytic_laplace.py` | transfer functions of the cascade; exact periodic scores, ramp and pulse rules, and their comparison with the numerical grid |
| `analytic/A_*.csv`–`E_*.csv`, `figS_analytic_vs_numeric.*` | exact versus numerical boundaries, convergence with time step, design rules; Figure S1 of the Supporting Information |
| `canonical_dense_boundary_grid_tau_sensor_1min.csv` | 71 transport times at a 1 min sensor constant; the grid plotted in Fig. 2c |
| `07_empirical_stress_test_full_1295.csv` | seven digitized traces against 185 parameter settings |
| `round2_analyses/A_*.csv` | envelopes at sensor constants of 0 and 1 min, binding criterion per cell, delay-aligned scores |
| `round2_analyses/B_*.csv` | coupled-volume budget and placement of a published 72 nL well |
| `round2_analyses/C_*.csv` | time-varying flux: steady, burst and declining at equal mean |
| `round2_analyses/D_*.csv` | fill--measure--clear and chronological-window results |
| `10_digitization_manifest.csv`, `11_digitization_qc.csv` | source locators, axis calibration and quality control for every digitized trace |

## Provenance of the empirical traces

The waveform proxies were digitized from figures in published articles. Each
trace retains its source identifier, figure and panel locator, axis calibration
points and digitization uncertainty. They are device outputs, not physiological
ground truth, and are redistributed here only as derived coordinates. Please
cite the original articles listed in `10_digitization_manifest.csv` when using
them.

## License

This repository carries two licenses.

| Path | License | File |
|---|---|---|
| `code/` | MIT | `LICENSE` |
| `outputs/` | CC BY 4.0 | `LICENSE-CC-BY-4.0.txt` |

You may reuse the code freely under the MIT terms, and share or adapt the data
files under CC BY 4.0 provided you give appropriate credit. Each directory also
carries a short `LICENSE-NOTICE.md` stating which license applies to it.

The empirical waveform files are coordinates digitized from figures in
published articles. They are redistributed here as derived coordinates with
full source locators; the underlying articles remain the copyright of their
publishers.

## Citing this archive

This archive is deposited at https://doi.org/10.5281/zenodo.22932238 (this DOI always resolves to the latest version). `CITATION.cff` carries the
machine-readable citation and `.zenodo.json` the deposit metadata. Please cite
both this archive and the Perspective it accompanies.

> Lu, Z.; Montagut Marques, M. J.; Hirose, K.; Umezu, S. Temporal fidelity in
> wearable sweat sensing: code and data, v1.1.0. Zenodo, 2026.
> https://doi.org/10.5281/zenodo.22932238

## Version history

- v1.1.0 (2026-11-02): Laplace-domain solutions of the transport–sensor cascade (`run_analytic_laplace.py`, `outputs/revision_2026-10-27/analytic/`); legends on every column of Figure 3c.
- v1.0.0 (2026-09-24): first release.
