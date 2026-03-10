# Kuramoto Model — Synchronization Dynamics

**Course:** Modelling Complex Systems, Module 2  
**Institute:** IISER Mohali  
**Author:** Avneesh Singh (MS23249)

## Overview

Numerical simulation of the all-to-all Kuramoto model with sine coupling.
The synchronization transition is characterized using the complex order
parameter r(t) for both Lorentzian and Gaussian frequency distributions.

## Model

The equation of motion for each oscillator:

$$\dot{\theta}_i = \omega_i + \frac{K}{N} \sum_{j=1}^{N} \sin(\theta_j - \theta_i)$$

Critical coupling from mean-field theory: $K_c = 2 / [\pi g(0)]$

## Results

![Phase snapshots](figures/new_fig1_phases.png)
![Time evolution](figures/new_fig2_timeseries.png)
![Synchronization transition](figures/new_fig3_transition.png)
![Finite-size effects](figures/new_fig4_fsize.png)

## How to Run
```bash
pip install numpy matplotlib
python code/kuramoto_sim.py
```

Figures are saved to the `figures/` directory.

## Report

Full report with theory and derivations: [`report/MS23249_m2_mcs.pdf`](report/kuramoto.pdf)
