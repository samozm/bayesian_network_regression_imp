# Simulations and other relevant materials for Ozminkowski and Solis-Lemus (2022)

Manuscript:
```
@article{Ozminkowski2022,
author = {Ozminkowski, S. and Sol'{i}s-Lemus, C.},
year = {2022},
title = {{Identifying microbial drivers in biological phenotypes with a Bayesian Network Regression model}},
journal = {In preparation}
}
```

## Folder structure:
- `simulation` folder contains the scripts to replicate the simulations in Ozminkowski and Solis-Lemus (2022).
- `additional_materials` folder contains the scripts used in the [JuliaCon2021 presentation](https://www.youtube.com/watch?v=ZYtyD8-Cweg) and the comparison to the R implementation in [Guha and Rodriguez (2018)](https://www.tandfonline.com/doi/full/10.1080/01621459.2020.1772079).


## Simulation

### Downloading the Julia package

To replicate the simulations, you need to download the Julia package [BayesianNetworkRegression.jl](https://github.com/samozm/BayesianNetworkRegression.jl). You 


1. Git clone the repository in your machine
```
git clone https://github.com/samozm/BayesianNetworkRegression.jl.git
```

2. Open julia and enter package mode by pressing `]`

3. Type the following command: 
```
dev PATH/BayesianNetworkRegression.jl
```
where PATH is the path to the `BayesianNetworkRegression.jl` repository in your machine.

### Generating the simulating data

All files for generating simulations are in the `simulation/` directory, but they need to be run outside the `simulation` folder.

Simulations are split into two different types: 'realistic' and 'unrealistic' (denoted 'theoretical' in the manuscript).

Unrealistic simulations are generated by running the following (from the main directory):
```
simulation/generate_unrealistic_sims.sh [SAMPLESIZE]
```

and realistic simulations are generated by running the following (from the main directory):
```
simulation/generate_realistic_sims.sh [SAMPLESIZE]
```
where \[SAMPLESIZE\] is the number of samples to generate.

### Understanding the output files

All the output files will be saved in two main folder:
- `data/simulation/unrealistic`
- `data/simulation/realistic`

For each combination of mu, pi, number of microbes (k), simulation type (unrealistic/realistic), and sample size, there are five generated files with the following file names:

- File `R=9_edge_mu=0.4_mu=[MU]_n_microbes=[K]_nu=10_out=bs_pi=[PI]_samplesize=[SAMPLESIZE]_simnum=2_type=[SIMTYPE].csv` contains the true B matrix (column matrix) associated with the simulated data.
- File `R=9_edge_mu=0.4_mu=[MU]_n_microbes=[K]_nu=10_out=main_effects_pi=[PI]_samplesize=[SAMPLESIZE]_simnum=2_type=[SIMTYPE].csv` contains the main effect for each node used to generate the results (since realistic simulations are generated using node main effects)
- File `R=9_edge_mu=0.4_mu=[MU]_n_microbes=[K]_nu=10_out=ms_pi=[PI]_samplesize=[SAMPLESIZE]_simnum=2_type=[SIMTYPE].csv` contains the m value (which shows which taxa were "sampled" for each "sample" (row in the X matrix/entry in the y vector)). A 1 indicates that taxa was included in that sample, a 0 indicates that it wasn't.
- File `R=9_edge_mu=0.4_mu=[MU]_n_microbes=[K]_nu=10_out=xis_pi=[PI]_samplesize=[SAMPLESIZE]_simnum=2_type=[SIMTYPE].csv` contais a vector of boolean values indicating whether each node is "turned on" - an edge has nonzero value if it connects two "turned on" nodes.
- File `R=9_edge_mu=0.4_mu=[MU]_n_microbes=[K]_nu=10_out=XYs_pi=[PI]_samplesize=[SAMPLESIZE]_simnum=2_type=[SIMTYPE].csv` contains the X matrices and corresponding y vectors for the simulation.


### Fitting the BNR model on simulated data

All scripts should be run at the main level `bayesian_network_regression_imp/`, not inside the `simulation` folder, and will assume that the output files have been stored in the two folders:
- `data/simulation/unrealistic`
- `data/simulation/realistic`

To fit the BNR model on the simulated data, we run 
```
julia --optimize=0 --math-mode=ieee --check-bounds=yes simulation/run_simulation_unrealistic.jl
```
for the unrealistic (theoretical) simulations, and
```
julia --optimize=0 --math-mode=ieee --check-bounds=yes simulation/run_simulation_realistic.jl
```
for the realistic simulations.


### Manuscript plots

The final plots in the manuscript can be obtained using the [`final-plots.Rmd`](https://github.com/samozm/bayesian_network_regression_imp/blob/main/additional_materials/final-plots.Rmd) file. Plots for the real data can be obtained using [`wagg_analysis/plots.Rmd](https://github.com/samozm/bayesian_network_regression_imp/tree/main/wagg_analysis/plots.Rmd) file. 


### Real Data Preparation and Augmentation

The data used for the real data section is in [`wagg_analysis/wagg_etal_data`](https://github.com/samozm/bayesian_network_regression_imp/tree/main/wagg_analysis/wagg_etal_data). Data preparation and augmentation can be found in [`prep_data.R`](https://github.com/samozm/bayesian_network_regression_imp/blob/main/wagg_analysis/src/prep_data.R).
