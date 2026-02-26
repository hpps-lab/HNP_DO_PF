# Hybrid Nonlinear Preconditioning for Double-Obstacle Phase-Field Simulations

This repository contains the source code used in the study:

S. Sakane, T. Takaki, "Hybrid nonlinear preconditioning approach for phase-field model with double-obstacle potential", Computational Materials Science 267 (2026) 114600. DOI: https://doi.org/10.1016/j.commatsci.2026.114600

![Graphical abstract](https://ars.els-cdn.com/content/image/1-s2.0-S0927025626001199-ga1_lrg.jpg)

This code implements a hybrid nonlinear preconditioning method that improves the numerical stability and accuracy of phase-field simulations using the double-obstacle (DO) potential, particularly under coarse interface resolution. The method combines nonlinear preconditioning models based on both double-obstacle and double-well (DW) potentials.

## Features

- Hybrid formulation of nonlinear preconditioning for the DO phase-field model.
- Supports stable and accurate simulations even with coarse spatial resolution.
- Demonstrated on:
  - Single grain growth
  - Triple junction migration under anisotropic grain boundary energies

## Repository Structure
```
HNP_DO_PF/
├── 2D_MPF_DO/                   # Conventional Implementation
│   ├── SingleGrain/                 # Sample for single grain growth (δ=4Δx, Δf=5e5 J/m^3)
│   │   ├── 2D_MPF_DO.c              # Source code 
│   │   ├── make.sh
│   │   └── run.sh
│   └── TripleJunction/              # Sample for triple junction migration (δ=4Δx, γA/γB=1.0)
│       ├── 2D_MPF_DO.c              # Source code
│       ├── make.sh
│       └── run.sh
├── 2D_MPF_HNP_DO/               # hybrid nonlinear preconditioning Implementation
│   ├── SingleGrain/                 # Sample for single grain growth (δ=4Δx, Δf=5e5 J/m^3)
│   │   ├── 2D_MPF_HNP_DO.c          # Source code
│   │   ├── make.sh
│   │   └── run.sh
│   └── TripleJunction/              # Sample for triple junction migration (δ=4Δx, γA/γB=1.0)
│       ├── 2D_MPF_HNP_DO.c          # Source code
│       ├── make.sh
│       └── run.sh
├── LICENSE
└── README.md                    # This file
```
## Requirements

- gcc

## Usage

After navigating to the each directory containing the source code, run the following command:
```bash
./make.sh
./run.sh
```
## Output Files

After running the simulation, the following output files are generated:
- ```computational_condition.dat``` file
  - Contents: computational_condition data in the simulation.
- ```.vtk``` files
  - Purpose: Store spatial distributions of phase-field variables at each time step.
  - Visualization: Can be visualized using ParaView or other VTK-compatible tools.
  - Location: Saved in the current directory, with one file per time step.

- ```interface_velocity_x_at_y00h.dat``` file
  - Contents: Time-series data recorded during the simulation.
  - Columns:
    1. Time step       [Δt]
    2. Simulation time [s]
    3. Interface position of grain 0 at y = 0 [m]
    4. Interface velocity of grain 0 at y = 0 [m/s]

- ```interface_velocity_x_at_y05h.dat``` file
  - Contents: Time-series data recorded during the simulation.
  - Columns:
    1. Time step       [Δt]
    2. Simulation time [s]
    3. Interface position of grain 0 at y = 0.5H [m]
    4. Interface velocity of grain 0 at y = 0.5H [m/s]


## License

This code is released under the MIT License. See LICENSE for details.


## Citation

If you use this approach in your research, please cite the following paper:

S. Sakane, T. Takaki, "Hybrid nonlinear preconditioning approach for phase-field model with double-obstacle potential", Computational Materials Science 267 (2026) 114600. DOI: https://doi.org/10.1016/j.commatsci.2026.114600
