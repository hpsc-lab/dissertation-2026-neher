# Reproducibility Repository

This repository contains all code to reproduce the numerical experiments and simulations presented in the dissertation. It includes validation scripts for pressure models and SPH methods, as well as patient-specific aortic simulations.

> **Note:** Data (`out/` and `data/` directory) will be archived on Zenodo after the dissertation is published. TODO

**Table of Contents**
- [System Requirements](#system-requirements)
- [Installation](#installation)
- [Repository Structure](#repository-structure)
- [Workflow Overview](#workflow-overview)
- [Running Simulations](#running-simulations)
- [Generating Figures](#generating-figures)

---

## System Requirements

### Required Software
- **Julia**: v1.12.6 (tested version) — [Install Julia](https://julialang.org/downloads/)
  - Other v1.12.x versions should work, but v1.12.6 guarantees reproducibility
- **ParaView**: 5.10+ (optional, only for interactive 3D visualizations)

### Hardware Requirements (recommended)

> **Note:** These are rough estimates based on typical runs. Your actual requirements may vary significantly depending on particle spacing, number of subjects, and HPC system configuration.

| Scenario | CPU Cores | RAM | Disk |
|----------|-----|-----|------|
| **Validation scripts only** | Any | ≥ 8 GB | ≥ 2 GB |
| **Patient simulations (CPU)** | ≥ 16 | ≥ 32 GB | ~50–200 GB |
| **Patient simulations (GPU)** | ≥ 8 | ≥ 16 GB | ~50–200 GB |

### GPU Support (Optional)

GPU scripts are available in [`scripts/aorta/gpu/`](scripts/aorta/gpu/) for faster execution.

**Backend Options:**
- `CUDABackend()` — NVIDIA GPUs
- `AMDGPUBackend()` — AMD GPUs
- `MetalBackend()` — Apple Silicon Macs

You **must adjust the backend** in GPU scripts to match your hardware. See [TrixiParticles.jl GPU documentation](https://trixi-framework.org/TrixiParticles.jl/stable/gpu/) for details.


---

## Installation

### 1. Install Julia

Install Julia following the instructions at https://julialang.org/downloads/.
Verify installation:
```bash
julia --version  # Should show v1.12.6+
```

### 2. Set Up Julia Environment
To install all necessary Julia packages, execute the following statement from within the folder that
contains the `README.md` file you are currently reading:
```bash
# Install all required packages
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```
This will recreate the exact Julia environment we used to obtain our results for full
reproducibility.


The project uses a custom package `SimulationSetup.jl` (located in [SimulationSetup/](SimulationSetup/)) for:
- Configuration management
- Logging and output formatting
- Preprocessing utilities
- Postprocessing helpers

---

## Repository Structure

```
code/
├── README.md                    # This file
├── Project.toml                 # Julia project configuration
├── Manifest.toml                # Locked dependencies (reproducibility)
├── SimulationSetup/             # Custom Julia package for utilities
│   ├── src/
│   │   ├── config.jl           # Configuration handling
│   │   ├── io.jl               # Input/output utilities
│   │   ├── logging.jl          # Logging setup
│   │   └── ...
│   └── Project.toml
├── scripts/
│   ├── validation/              # Validation tests for models
│   │   ├── pressure_model/      # Windkessel & vessel models
│   │   ├── open_boundaries/     # SPH boundary condition tests
│   │   └── fsi/                 # Fluid-structure interaction tests
│   ├── aorta/                   # Patient-specific aortic simulations
│   │   ├── setup_*.jl           # Setup scripts
│   │   ├── simulate_*.jl        # Simulation drivers
│   │   └── gpu/                 # GPU-accelerated versions
│   ├── pbs_job_scripts/         # HPC job submission scripts
│   ├── visualization/           # Figure generation scripts
│   └── ...
└── data/
    ├── aorta_centered/          # Segmented aorta geometries
    ├── aorta_preprocessed/      # Pre-processed geometries
    ├── reference_data/          # Validation reference data
    └── ...
```

---

## Workflow Overview

### **Chapter 2: Pressure Models**
- Windkessel pressure models (RC, RCR)
- Validation against experimental data
- **Time:** Minutes

### **Chapter 3: SPH Method**
- Smoothed Particle Hydrodynamics implementation
- Open boundary conditions
- Fluid-structure interaction (FSI) coupling
- **Time:** About an hour

### **Chapter 4: Preprocessing Tool**
Independent reproducibility repository at [Zenodo](https://zenodo.org/records/17384814)

### **Chapter 5: Patient-Specific Aortic Simulations** (*Computationally intensive*)
- 14 patient subjects (F01–F16, excluding F04, F06)
- 3 hemodynamic scenarios (normotensive, exercise, hypertensive)
- 2 models (rigid, elastic)
- Requires **preprocessing → warmup → production** workflow
- **Time:** Hours to days per subject (depending on hardware)

---

## Running Simulations

### Chapter 2: Pressure Models


**Prerequisites:** None (data included)

```bash
# Validate Windkessel pressure models (RC and RCR)
julia --project=. scripts/validation/pressure_model/windkessel_model.jl

# Model pressure dynamics in healthy, stenosed, and arteriosclerotic vessels
julia --project=. scripts/validation/pressure_model/vessel_model.jl
```

---

### Chapter 3: SPH Method

**Prerequisites:** None (data included)

```bash
# Poiseuille flow 2D (steady-state reference)
julia --project=. scripts/validation/open_boundaries/poiseuille_flow_2d.jl

# Pulsatile channel flow 3D (time-varying flow)
julia --project=. scripts/validation/open_boundaries/validation_pulsatile_channel_flow_3d.jl

# Pulse wave propagation with FSI coupling
julia --project=. scripts/validation/fsi/pulse_wave_propagation_3d.jl
```

**Output:** Results saved to `out/validation/`

---

### Chapter 4: Preprocessing Tool

Full reproducibility repository: https://zenodo.org/records/17384814

---

### Chapter 5: Patient-Specific Aortic Simulations

**This section requires:**
- HPC system with PBS job scheduler
- 50+ GB disk space
- Patience: Simulations take multiple hours/days depending on hardware

**Available Subjects:** F01–F16 (excluding F04, F06)

**Hemodynamic Scenarios:** normotensive, exercise, hypertensive

**Models:** rigid (faster) or elastic/FSI (more realistic, slower)

#### **Workflow Dependency**

All three stages must be run sequentially for each subject:

```
Stage 1: Preprocessing
    ↓ (Creates: data/aorta_initial_condition/{subject}/)
Stage 2: Warmup Simulation (7 cardiac cycles)
    ↓ (Creates: out/aorta/{subject}/warmup/)
Stage 3: Production Simulation (1 cardiac cycle)
    ↓ (Creates: out/aorta/{subject}/cycle/)
```

**Do not skip Stage 2** — it establishes the periodic steady state required for Stage 3.

#### Stage 1: Preprocessing

Loads STL segmentations and generates simulation-ready particle distributions.

**Parameters:**
- `SUBJECT`: Patient ID (default: F09) — available: F01-F16 (excluding F04, F06)
- `PARTICLE_SPACING`: SPH particle spacing in meters (default: 0.001)
- `JULIA_THREADS`: Number of Julia threads (default: 128)
- `VERSION`: Output version label (default: v1.0.1)

```bash
# Default subject (F09)
qsub scripts/pbs_job_scripts/preprocessing.pbs

# Custom subject with finer spacing
qsub -v SUBJECT=F10,PARTICLE_SPACING=0.0005 \
     scripts/pbs_job_scripts/preprocessing.pbs
```

**Expected Output:**
```
data/aorta_initial_condition/v1.0/packed_results_F09/
```

> **Optional Step:** Run `scripts/aorta/preprocess_geometries.jl` if you want to process (centering and scaling) new/custom STL geometries. Pre-processed geometries for all subjects (F01–F16) are already included in `data/aorta_preprocessed/`.

---

#### Stage 2: Warmup Simulation

Runs 7 cardiac cycles to establish periodic steady state (required for production run).

> **This takes multiple hours per subject** depending on model and hardware. Use `RESTART=true` to resume from checkpoints if interrupted.

**Parameters:**
- `MODEL`: `rigid` or `elastic`
- `SCENARIO`: `normotensive`, `exercise`, or `hypertensive`
- `DEVICE`: `cpu` or `gpu`
- `RESTART`: `true` to resume from checkpoint, `false` to start fresh
- `SUBJECT`: Patient ID (default: F10)
- `PARTICLE_SPACING`: Must match preprocessing! (default: 0.001)
- `JULIA_THREADS`: Number of threads (default: 128)
- `VERSION`: Version label (default: v1.0.1)

```bash
# Default: Rigid model on CPU
qsub scripts/pbs_job_scripts/transient_aorta.pbs

# FSI model on GPU, hypertensive scenario
qsub -v MODEL=elastic,DEVICE=gpu,SCENARIO=hypertensive,SUBJECT=F09 \
     scripts/pbs_job_scripts/transient_aorta.pbs

# Resume from checkpoint
qsub -v RESTART=true,MODEL=elastic,SUBJECT=F10 \
     scripts/pbs_job_scripts/transient_aorta.pbs
```

**Expected Output:**
```
out/out_normotensive/F10/rigid/
```

---

#### Stage 3: Production Simulation

Runs 1 cardiac cycle from the periodic state established in Stage 2.

**Prerequisites:** Stage 2 must have completed for your subject/model/scenario combination

**Parameters:** Same as Stage 2 (must match!)

```bash
# Rigid model on CPU (default)
qsub scripts/pbs_job_scripts/simulate_cycle_aorta.pbs

# Elastic FSI model on GPU, exercise scenario
qsub -v MODEL=elastic,DEVICE=gpu,SCENARIO=exercise,SUBJECT=F09 \
     scripts/pbs_job_scripts/simulate_cycle_aorta.pbs

# Resume from checkpoint
qsub -v RESTART=true,MODEL=elastic,SUBJECT=F10 \
     scripts/pbs_job_scripts/simulate_cycle_aorta.pbs
```


**Expected Output:**
```
out/out_normotensive/F10/rigid/full_cycle/
```

---

## Generating Figures

After running simulations, generate figures using visualization scripts.

**Two options:**
1. **Use pre-computed results:** Pre-processed results are in `data/` and `out/`
2. **Use your own results:** Re-run simulations first, then point visualization scripts to your output

### Chapter 2: Pressure Models

```bash
# Fig. 2.6 (a) and (b) — Windkessel model validation
julia --project=. scripts/visualization/pressure_model/windkessel_model.jl
# → Output: figures/windkessel_model/

# Fig. 2.7 — Vessel model under pulsatile flow
julia --project=. scripts/visualization/pressure_model/vessel_model.jl
# → Output: figures/windkessel_model/
```

---

### Chapter 3: SPH Method

```bash
# Kernel functions
# Fig. 3.4
julia --project=. scripts/visualization/sph/kernel_2d.jl
# → Output: figures/sph/

# Fig. 3.5
julia --project=. scripts/visualization/sph/kernel_1d.jl
# → Output: figures/sph/
```

```bash
# Fig. 3.10 — Mirroring methods comparison
julia --project=. scripts/visualization/open_boundaries/mirroring_methods.jl
# → Output: figures/open_boundaries/

# Fig. 3.16 — Kernel ramping technique
julia --project=. scripts/visualization/open_boundaries/kernel_ramping.jl
# → Output: figures/open_boundaries/
```

```bash
# Fig. 3.20 — Poiseuille flow 2D
julia --project=. scripts/visualization/validation/poiseuille_flow_2d.jl
# → Output: figures/validation/

# Fig. 3.21 (a) — Wall shear stress (WSS) analysis, time point 0
julia --project=. scripts/visualization/validation/wss_poiseuille_flow_2d.jl 0
# → Output: figures/validation/

# Fig. 3.21 (b) — WSS analysis, time point 1
julia --project=. scripts/visualization/validation/wss_poiseuille_flow_2d.jl 1
# → Output: figures/validation/

# Fig. 3.22 — Pulsatile channel flow 3D
julia --project=. scripts/visualization/validation/pulsatile_channel_flow_3d.jl
# → Output: figures/validation/

# Fig. 3.26 — Pulse wave propagation with FSI
julia --project=. scripts/visualization/validation/pulse_wave_propagation_3d.jl
# → Output: figures/validation/
```

**Interactive 3D visualizations (ParaView):**

```bash
# Fig. 3.15 — Robustness test results
paraview scripts/visualization/paraview_states/robustness.pvsm

# Fig. 3.24 — Pulsatile pipe setup visualization
paraview scripts/visualization/paraview_states/setup_pulsatile_pipe.pvsm

# Fig. 3.25 — Pulsatile pipe 32-particle result
paraview scripts/visualization/paraview_states/pulsatile_pipe_32.pvsm
```

---

### Chapter 5: Patient-Specific Aortic Simulations

```bash
# Fig. 5.9 (a) and 5.10 — Windkessel convergence across particle spacings
julia --project=. scripts/visualization/aorta/wk_convergence.jl
# → Output: figures/aorta/

# Fig. 5.9 (b) — Normalized convergence rates
julia --project=. scripts/visualization/aorta/wk_convergence_normalized.jl
# → Output: figures/aorta/

# Fig. 5.11 (a) and (b) — Periodic state under hypertensive scenario
julia --project=. scripts/visualization/aorta/wk_periodic_state.jl F09 1 hypertensive
julia --project=. scripts/visualization/aorta/wk_periodic_state.jl F10 1 hypertensive
# → Output: figures/aorta/

# Fig. 5.11 — Periodic state under exercise scenario
julia --project=. scripts/visualization/aorta/wk_periodic_state.jl F09 1 exercise
# → Output: figures/aorta/

# Fig. 5.13 (a) and (b) — Windkessel model comparison across subjects
julia --project=. scripts/visualization/aorta/wk_comparison.jl F09 0
julia --project=. scripts/visualization/aorta/wk_comparison.jl F10 0
# → Output: figures/aorta/

# Fig. 5.14 — FSI vs. rigid model comparison
julia --project=. scripts/visualization/aorta/wk_comparison_fsi.jl F09
# → Output: figures/aorta/

# Fig. 5.15 — Volume flow rate comparison
julia --project=. scripts/visualization/aorta/volume_rate_comparison.jl F09
# → Output: figures/aorta/

# Fig. 5.16 — Branch-wise flow distribution
julia --project=. scripts/visualization/aorta/branchwise_flow_comparison.jl F09
# → Output: figures/aorta/
```
