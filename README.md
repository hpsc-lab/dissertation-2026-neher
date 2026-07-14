# A Particle-Based Simulation Pipeline for Patient-Specific Aortic Hemodynamics

This repository contains information and code to reproduce the results presented in the
dissertation
```bibtex
TODO
@misc{neher2026diss,
  author    = {Neher, Niklas S.},
  title     = {A Particle-Based Simulation Pipeline for Patient-Specific Aortic Hemodynamics}
}
```

If you find these results useful, please cite the work mentioned above. If you use the implementations provided here, please also cite this repository as
```bibtex
@software{neher2026reproducibility,
  author       = {Neher, Niklas S.},
  title        = {Reproducibility repository for
                  "A Particle-Based Simulation Pipeline for Patient-Specific Aortic Hemodynamics"},
  month        = jul,
  year         = 2026,
  publisher    = {Zenodo},
  version      = {v0.1},
  doi          = {10.5281/zenodo.21357207},
  url          = {https://doi.org/10.5281/zenodo.21357207},
}
```

> Note: This thesis has been submitted for a doctoral degree at the Faculty of Mathematics, Natural Sciences, and Engineering of the University of Augsburg. It is currently under review  and has not yet been formally accepted or published.

## Abstract

This thesis develops a particle-based simulation pipeline for patient-specific hemodynamic analysis of the human aorta.
Since detailed *in vivo* measurements of clinically relevant hemodynamic quantities are difficult and often invasive,
*in silico* simulations provide an important complement to morphology-based assessment of vascular disease.
However, such simulations present significant engineering challenges, particularly regarding the complexity of anatomical geometries, appropriate boundary conditions, and the robust coupling of flow and structure.
To address these challenges, this work adopts a particle-based approach employing Smoothed Particle Hydrodynamics (SPH) and presents a comprehensive, largely automated simulation framework.

The work introduces two main contributions.
First, a fully particle-based preprocessing methodology is developed that automatically converts complex geometries into simulation-ready initial conditions.
This technique utilizes a face-based neighborhood search to construct a memory-efficient Signed Distance Field (SDF), enabling localized computations near surface regions.
To create an initial particle configuration, a hierarchical winding number method for fast and accurate inside-outside segmentation is applied.
Particle positions are then relaxed using an SPH-inspired scheme, where the SDF acts as a geometric constraint to prevent particles from drifting outside the domain.
By incorporating dynamic packing of boundary particles, this approach ensures an accurate and robust representation of complex geometry surfaces through purely local particle interactions, fully leveraging the meshless nature of the method.
Second, building upon this preprocessing tool, a generalized simulation pipeline for patient-specific aortic hemodynamics is established within the TrixiParticles.jl framework.
The pipeline integrates geometry processing, numerical modeling, and Fluid Structure Interaction (FSI).
To evaluate the practical applicability of the framework, simulations were conducted on different patient-specific geometries under various physiological scenarios, comparing rigid and compliant vessel wall models.

The results demonstrate that the developed framework yields physiologically plausible hemodynamic behavior.
In particular, compliant wall models show more realistic flow dynamics, as reflected in Wall Shear STress (WSS) analyses.
Overall, this work provides a robust, modular, and highly automated computational foundation that bridges the gap between purely morphological assessments and quantitative hemodynamic predictions.
By design, the framework lays the groundwork not only for future clinical validation but also for the seamless integration of upcoming physical and biological model extensions.


## Numerical experiments

The numerical experiments presented in the dissertation use
[TrixiParticles.jl](https://github.com/trixi-framework/TrixiParticles.jl).
To reproduce the numerical experiments, you need to install
[Julia](https://julialang.org/).

The subfolder [`code`](/code/) of this repository contains a `README.md` file with
instructions to reproduce the numerical experiments.
The subfolders also include the input data, result data and scripts for postprocessing.

All numerical experiments were carried out using Julia v1.12.6.

## License

The contents of this repository are available under the [MIT license](LICENSE.md). If you reuse our
code or data, please also cite me (see above).

## Disclaimer

Everything is provided as is and without warranty. Use at your own risk!
