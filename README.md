# Does the 95% herd immunity threshold for measles hold under spatial heterogeneity?

**A two-patch SEIR-V model for measles, with application to South African district data**

Jemma Antoncich (26927748) and Michael Zietsman (27210316)
Biomathematics 374, Department of Mathematical Sciences, Stellenbosch University, 2026

[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=jemmaantoncich-web/Measles-herd-immunity&file=measles_model.m)

---

## What this code does

`measles_model.m` reproduces every number, table and figure in the Results section and Appendix of our report. It:

1. solves the two-patch SEIR-V model (equations 1–8 of the report) for three vaccination scenarios;
2. calculates the control reproduction number R_c with the next generation matrix (equations 11–12);
3. calculates R_c across all mixing fractions 0 ≤ ϕ ≤ 1;
4. calculates the normalised sensitivity indices of R_c (equation 13);
5. finds the second-dose coverage in Capricorn needed for R_c = 1;
6. repeats the simulations with R₀ = 20.

## How to run it

**Option 1: MATLAB Online (no installation).**
Click the **Open in MATLAB Online** button above, sign in with a MathWorks account, and type in the Command Window:

```matlab
measles_model
```

**Option 2: MATLAB on your computer.**
Download this repository (green **Code** button, then **Download ZIP**), unzip it, open the folder in MATLAB, and type:

```matlab
measles_model
```

The script runs in under a minute. No input or setup is needed.

### Requirements

- MATLAB **R2020a or later** (the code uses `exportgraphics`, `yline` and `table`).
- **No toolboxes** are needed. Only built-in MATLAB functions are used (`ode45`, `fzero`, `trapz`, `eig`).

## Output

All output files are saved in the current folder.

| Output | In the report |
|---|---|
| Command Window: Table 5 summary | Table 5 |
| Command Window: R₀ = 20 check | Section 4.2 |
| `Figure2.png` – infectious proportion, all scenarios | Figure 2 |
| `Figure3.png` – susceptible and vaccinated proportions | Figure 3 |
| Command Window: infections by vaccination status | Table 6 |
| `Figure4.png` – R_c against mixing fraction ϕ | Figure 4 |
| `Figure5.png` – sensitivity indices | Figure 5 |
| Command Window: critical coverage | Section 4.5 |
| `FigureA1.png` to `FigureA3.png` – all compartments per scenario | Figures A1–A3 |
| `FigureA4.png` – infectious proportion with R₀ = 20 | Figure A4 |
| `results.txt` – copy of everything printed in the Command Window | All tables |

## Structure of the code

All functions are in the single file `measles_model.m`.

| Function | Purpose |
|---|---|
| `measles_model` | Main script: sets parameters and scenarios, runs every analysis, prints tables and saves figures |
| `seir_v_odes` | The system of ODEs, equations (1)–(8) |
| `calc_Rc` | Builds the matrices F and V and returns R_c as the spectral radius of FV⁻¹ (equation 11) |
| `calc_R0` | Basic reproduction number R₀ = βσ / ((σ + μ)(γ + μ)) |
| `run_simulation` | Sets the initial conditions (one exposed person in Patch 1) and solves the ODEs with `ode45` for 730 days |
| `outbreak_size` | Peak infectious, peak day, day the outbreak ends and total infected (cumulative new infections) |
| `infections_by_status` | Splits total infections into unvaccinated, one-dose and two-dose individuals |
| `sensitivity` | Normalised forward sensitivity indices (equation 13), by central difference |
| `plot_infectious`, `plot_all_compartments` | Plotting functions for Figure 2, Figure A4 and Figures A1–A3 |

## Parameters and scenarios

| Parameter | Value | Meaning |
|---|---|---|
| β | 1.25 day⁻¹ | Transmission rate (gives R₀ ≈ 10) |
| σ | 0.125 day⁻¹ | Progression rate E → I |
| γ | 0.125 day⁻¹ | Recovery rate I → R |
| μ | 0.0000418 day⁻¹ | Birth/death rate, 1/(65 × 365) |
| ϵ₁, ϵ₂ | 0.85, 0.97 | Efficacy of one and two doses |
| ϕ | 0.10 | Mixing fraction between patches |

| Scenario | Patch 1 (v₁, v₂) | Patch 2 (v₁, v₂) | Populations |
|---|---|---|---|
| 1: Uniform | 0.05, 0.90 | 0.05, 0.90 | 400 000 and 1 600 000 |
| 2: Clustered | 0.20, 0.60 | 0.0475, 0.94 | 400 000 and 1 600 000 |
| 3: Capricorn–Johannesburg | 0.081, 0.755 | 0.045, 0.904 | 1 447 103 and 4 803 262 |

To change a parameter, edit Section 1 (base parameters) or Section 2 (scenarios) at the top of `measles_model.m` and run the script again.

## Data sources

- Vaccination coverage: National Department of Health (2020), *Expanded Programme on Immunisation (EPI) National Coverage Survey Report 2020* (2019 survey data).
- District populations: Statistics South Africa (2023), *Census 2022*.
- Full references are in the report.

## Generative AI use

Generative AI (Claude, Anthropic) was used to help write, restructure and check this code. All code was run and checked by the authors. See the generative AI use declaration in the report.

