# Social Determinants and the Identifiability Problem in Computational Epidemiology

Code, data, results, and figures accompanying the article **"Social Determinants and the Identifiability Problem in Computational Epidemiology."**

## Overview

Social determinants of health can change epidemic dynamics through different mechanisms. Two of them are increased exposure (higher transmission) and reduced access to care (slower recovery). This article shows that these mechanisms can be **practically unidentifiable** from early outbreak case counts, which is the data typically available at the start of an epidemic.

The analysis uses a SEIR-type compartmental model for SARS-CoV-2 with asymptomatic ($I_A$) and symptomatic ($I_S$) infectious classes, and has three parts:

1. **Analytical degeneracy.** For the early exponential growth rate $r$ there is a curve $\beta(\nu, r)$ of transmission/recovery pairs that all give exactly the same growth rate. The same holds for $\mathcal{R}_0$. Derivations are in the Supplementary Material (S3).
2. **Synthetic outbreak.** A known-mechanism outbreak is generated with Poisson observation noise. Two competing hypotheses are each fit to it: **High Transmission (A)** with $\nu = 0.30$ and **Low Recovery (B)** with $\nu = 0.05$. Both estimate $\beta$.
3. **Real outbreak.** The same procedure is applied to early US COVID-19 daily case counts (21 January – 30 March 2020).

In both cases the two hypotheses reproduce the data almost identically (ΔAIC is small), even though they imply different social mechanisms and different interventions.

## Repository structure

```
.
├── Data/
│   ├── PrepareOWIDDataV5.m                    # Filters raw OWID/ECDC data to the US series and the 70-day window
│   ├── full_raw_data.csv                      # Archived OWID/ECDC compilation (all locations)
│   ├── owid_archived_ecdc_us_daily_cases.csv  # US daily cases, 31 Dec 2019 – 30 Jun 2020
│   ├── real_data_window.csv                   # 70-day analysis window (21 Jan – 30 Mar 2020)
│   └── Synthetic_Observations.csv             # Synthetic outbreak: true incidence + Poisson-noised observations
├── Synthetic/
│   └── RunnerSyntheticV5.m                    # Synthetic experiment (Fig. 2, Table 1, 200-replicate check)
├── RealData/
│   └── RunnerRealDataV5.m                     # Real-data experiment (Fig. 3, Table 1, window sensitivity)
├── Supporting/
│   └── SupInforV5.m                           # Supplementary figures and the β_A ≠ β_S extension
├── Figures/                                   # Main-text and supplementary figures (PNG)
│   └── Sup/
└── Results/                                   # Numerical results (CSV) and LaTeX tables (.tex)
```

## Requirements

- MATLAB. `exportgraphics` requires **R2020a or later**.
- No additional toolboxes are called explicitly. The code uses base-MATLAB `ode45`, `fminbnd`, and `fminsearch`. `poissrnd`, used to generate synthetic noise, is part of the Statistics and Machine Learning Toolbox.
- The scripts locate the project root automatically by searching upward from the current folder for a `Data/` directory. Run them from anywhere inside the repository.

## How to reproduce

Run the scripts in this order:

| Step | Script | Purpose |
|---|---|---|
| 0 (optional) | `Data/PrepareOWIDDataV5.m` | Rebuilds `owid_archived_ecdc_us_daily_cases.csv` and `real_data_window.csv` from `full_raw_data.csv`. The filtered files are already included. |
| 1 | `Synthetic/RunnerSyntheticV5.m` | Loads `Data/Synthetic_Observations.csv`, or generates it (seed `2026`) if missing. Fits both hypotheses and runs 200 noise replicates (seed `12345`). |
| 2 | `RealData/RunnerRealDataV5.m` | Estimates $\hat r$, fits both hypotheses to the real data, and runs the window-sensitivity analysis. |
| 3 | `Supporting/SupInforV5.m` | Produces the supplementary degeneracy figures and the two-parameter ($\beta_A$, $\beta_S$) fits. Needs `Synthetic_Observations.csv` from step 1. |

Outputs are written to `Figures/` and `Results/`. Existing files are overwritten.

To regenerate the synthetic dataset from scratch, delete `Data/Synthetic_Observations.csv` before step 1.

## Model and fitting

**Model.** SEAIR-type ODE system with demographic turnover and disease-induced mortality (Supplementary Material, Section S1). Parameters and initial conditions follow Ogbunugafor et al. (2020) and Gomez et al. (2020), and are listed in Supplementary Material, Section S2.

**Growth-rate degeneracy** (with $\beta_A = \beta_S = \beta$):

$$
\beta(\nu, r) = \frac{(\mu+\varepsilon+r)(\mu+\omega+r)(\mu_S+\nu+r)}{\varepsilon\left[(\mu_S+\nu+r) + (1-p)\,\omega\right]}
$$

**Observation model.** Model predictions are daily new symptomatic cases, obtained from the cumulative symptomatic incidence $C_S(t)$ with $dC_S/dt = (1-p)\,\omega\,I_A$. There are $n = 70$ observations per fit.

**Estimation.**
- One-parameter fits ($\beta$ only, $\nu$ fixed per hypothesis) minimize the sum of squared errors with `fminbnd` on $0.01 \le \beta \le 3$.
- The early growth rate $\hat r$ comes from a Poisson regression $y_t \sim \text{Poisson}(e^{a + rt})$ fitted by maximum likelihood (`fminsearch`). It starts from the day cumulative cases first reach 5.
- Two-parameter fits ($\beta_A$, $\beta_S$) use `fminsearch` from 12 starting points, with the bounds $0 \le \beta \le 5$ enforced by a sine reparameterization. The lowest-SSE solution is kept.
- Model comparison uses the Gaussian-likelihood AIC, $\text{AIC} = 2k - 2\log L$.

## Results

### Main analyses (shared β)

| Experiment | r (day⁻¹) | ν | β(ν, r) | Fitted β | R₀ | AIC |
|---|---|---|---|---|---|---|
| Synthetic truth | 0.1000 | 0.12 | 0.4897 | not estimated | 1.7043 | not estimated |
| Synthetic — High Transmission (A) | 0.1000 | 0.30 | 0.5032 | 0.5059 | 1.6514 | 498.85 |
| Synthetic — Low Recovery (B) | 0.1000 | 0.05 | 0.4765 | 0.4749 | 1.8861 | 500.48 |
| Real data — High Transmission (A) | 0.2079 | 0.30 | 0.7690 | 0.6360 | 2.0760 | 1160.43 |
| Real data — Low Recovery (B) | 0.2079 | 0.05 | 0.7495 | 0.6092 | 2.4191 | 1161.02 |

Source files: `Results/SyntheticSummaryTable.csv` and `Results/RealDataSummaryTable.csv`.

### Robustness checks

- **200 synthetic noise replicates** (`Results/SyntheticReplicates.csv`, `SyntheticReplicateSummary.csv`). The mean ΔAIC (A − B) is 0.60 with SD 4.32. High Transmission is favored in 93 of 200 replicates. |ΔAIC| is below 2 in 66 replicates, below 4 in 128, and above 10 in 5. The preferred hypothesis therefore depends on the noise realization.
- **Window sensitivity, real data** (`Results/RealDataWindowSensitivity.csv`). Shifting the 70-day window start by −2 to +4 days moves $\hat r$ between 0.17 and 0.24 per day and ΔAIC between −3.4 and +2.7. The sign of the preferred hypothesis changes, and every |ΔAIC| stays below 4.
- **β_A ≠ β_S extension** (`Results/betaAB_TwoParamResults.csv`). The two hypotheses differ by less than 0.1 AIC unit on both datasets.

| Dataset | Hypothesis | β̂_A | β̂_S | AIC |
|---|---|---|---|---|
| Synthetic | High Transmission | 0.505 | 0.533 | 500.84 |
| Synthetic | Low Recovery | 0.500 | 0.258 | 500.85 |
| COVID-19 | High Transmission | 0.662 | 0.000 | 1161.94 |
| COVID-19 | Low Recovery | 0.662 | 0.000 | 1161.94 |

For the real data, $\hat\beta_S$ reaches the lower bound of the search range under both hypotheses. At $\beta_S = 0$, $\nu$ drops out of the fitted incidence, so the two fits coincide exactly.

## Figures

| File | Description |
|---|---|
| `Figures/Figure2A_Synthetic_Degeneracy.png` | Degeneracy curve for r = 0.10, with the two hypotheses marked |
| `Figures/Figure2B_Synthetic_Fits.png` | Synthetic outbreak and the two fitted trajectories |
| `Figures/Figure3A_RealData_Degeneracy.png` | Degeneracy curve at the growth rate estimated from US data |
| `Figures/Figure3B_RealData_Fits.png` | Real-data fits under both hypotheses |
| `Figures/FigureS1_R_Degeneracy_Multi.png` | Growth-rate degeneracy curves for several values of r |
| `Figures/FigureS2_R0_Degeneracy_Multi.png` | $\mathcal{R}_0$ degeneracy curves |
| `Figures/FigureS2_RealData_SideBySide.png` | Independent real-data fits, shown side by side |
| `Figures/FigureS5_RealData_TwoParamFit.png` | Real data with independent β_A and β_S |
| `Figures/Sup/FigureS4_Synthetic_TwoParamFit.png` | Synthetic data with independent β_A and β_S |
| `Figures/Sup/SyntheticData.png` | The synthetic outbreak (true incidence and Poisson-noised observations) |
| `Figures/Sup/RealDataOnly.png` | The real case series |


## Data sources

Real case counts come from the archived ECDC-based compilation maintained by [*Our World in Data*](https://github.com/owid/covid-19-data). The series is daily newly confirmed US cases from 21 January to 30 March 2020, used without rescaling. `Data/full_raw_data.csv` is the archived file for all locations. The US series and the 70-day window are extracted from it by `PrepareOWIDDataV5.m`. The synthetic data are generated by `RunnerSyntheticV5.m` and contain no real-world data.
