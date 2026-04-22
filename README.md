# Computational Characterization of Inertial Particle Transport Regimes in Channel Flow

Source code for the senior thesis by Ian MacIver, Department of Mechanical and Aerospace Engineering, Princeton University, 2026.

Adviser: Michael E. Mueller | Second reader: Marcus Hultmark

## Overview

This repository contains the MATLAB code used to generate transport regime maps of inertial particles injected perpendicularly from one wall of a channel into a crossflow. The algorithm sweeps a two-dimensional parameter space defined by the Stokes number (St) and normalized injection velocity (v₀/U_ref), computing the maximum wall-normal penetration depth for each case. The resulting regime maps identify the boundaries between near-wall confinement, partial penetration, and midplane crossing.

The code operates in two phases:

- **Laminar phase:** Uses the analytical Poiseuille velocity profile as the carrier field. Validated against the exact penetration scaling y_max = v₀ · τ_p to within 0.14%.
- **Turbulent phase:** Replaces the Poiseuille profile with frozen snapshots of the JHTDB channel5200 DNS field at Re_τ = 5200. Eleven snapshots are swept independently and aggregated into an ensemble mean and variance.

A secondary gravity extension verifies the non-dimensional regime map translation analysis developed in the thesis.

## Requirements

- MATLAB R2024b (or later)
- Parallel Computing Toolbox (for `parfor` in the turbulent sweep)
- Internet connection (for JHTDB REST API queries during the turbulent phase)
- A JHTDB authentication key (free registration at [http://turbulence.pha.jhu.edu](http://turbulence.pha.jhu.edu))

## Repository Structure

```
├── README.md
├── LICENSE
├── laminar/
│   └── Final_Laminar_Flow_Model_Plot.m      % Laminar sweep + all laminar figures
├── turbulent/
│   └── run_turbulent_particle_sweep_Final.m % Single-snapshot turbulent sweep
├── postprocessing/
│   └── compute_ensemble_statistics.m        % Ensemble aggregation + turbulent figures
└── figures/                                 % Generated figures (created at runtime)
```

## Usage

### Laminar sweep

```matlab
cd laminar
Final_Laminar_Flow_Model_Plot
```

Runs the full 41 × 41 parameter sweep over laminar Poiseuille flow, generates all laminar figures (L1–L5), and prints the R5 analytical consistency verification. Runtime is approximately 1–2 minutes on a modern workstation. No external dependencies beyond base MATLAB.

### Turbulent sweep (single snapshot)

```matlab
cd turbulent
run_turbulent_particle_sweep_Final
```

Before running, open the script and set `turb.authkey` to your own JHTDB authentication key and set `turb.time0` to the desired snapshot (integer from 1 to 11). On first run, the script queries the JHTDB REST API for the velocity field on a 96 × 72 × 64 grid (442,368 points) and caches the result to disk as `channel5200_flow_cache_<k>.mat`. Subsequent runs with the same configuration load from cache. Each snapshot sweep completes in approximately 3–6 minutes using MATLAB's parallel pool.

To generate the full eleven-snapshot ensemble, run the script eleven times with `turb.time0` set to 1 through 11, producing eleven result files (`turbulent_particle_final_results_1.mat` through `_11.mat`).

### Ensemble aggregation

```matlab
cd postprocessing
compute_ensemble_statistics
```

Loads all eleven result files, validates parameter consistency, stacks the y_max maps into a [41 × 41 × 11] array, and produces the ensemble mean and standard deviation maps, per-snapshot boundary fits, convergence assessment, and all turbulent figures.

## Key Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| τ_p range | 10⁻³ to 10¹ outer units | 41 points, logarithmically spaced |
| v₀ range | 0.1 to 6.0 U_ref | 41 points, uniformly spaced |
| p_stick | 0.20 | Sticking probability at wall contact |
| N_bounce_max | 50 | Maximum elastic reflections before forced deposition |
| ε_wall | 10⁻⁶ outer units | Post-reflection wall offset |
| τ_p stiffness threshold | 0.1 outer units | ode15s below, ode45 above |
| RelTol / AbsTol | 10⁻⁶ / 10⁻⁹ | ODE solver tolerances (turbulent) |
| t_max | 50 outer units | Integration horizon per bounce segment |
| Grid resolution | 96 × 72 × 64 | JHTDB query grid (N_x × N_y × N_z) |

## JHTDB Authentication

The turbulent phase queries velocity data from the Johns Hopkins Turbulence Database channel5200 dataset. You will need a free authentication key:

1. Register at [http://turbulence.pha.jhu.edu](http://turbulence.pha.jhu.edu)
2. Copy your key from the account page
3. Replace the `turb.authkey` value in `run_turbulent_particle_sweep_Final.m`

## Citing This Work

If you use this code, please cite the thesis:

```
MacIver, I. (2026). Computational Characterization of Inertial Particle Transport
Regimes in Channel Flow. Senior thesis, Department of Mechanical and Aerospace
Engineering, Princeton University.
```

## License

MIT License. See [LICENSE](LICENSE) for details.
