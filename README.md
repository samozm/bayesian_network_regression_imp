#### Simulations and other relevant materials for "Identifying microbial drivers in biological phenotypes with a Bayesian Network Regression model".

All code is meant to be used with [BayesianNetworkRegression.jl](https://github.com/samozm/BayesianNetworkRegression.jl).

### Simulation
All files for generating simulations are in the `simulation/` directory.

Simulations are split into two different types: 'realistic' and 'unrealistic'.

Unrealistic simulations are generated by running the following (from the main directory):

```
simulation/generate_unrealistic_sims.sh [SAMPLESIZE]
```

and realistic simulations are generated by running the following (from the main directory):

```
simulation/generate_realistic_sims.sh [SAMPLESIZE]
```

where \[SAMPLESIZE\] is the number of samples to generate.

For each combination of μ,π,number of microbes, simulation type, and sample size five files are generated in the following format:


`R=9_edge_mu=0.4_mu=[MU]_n_microbes=[K]_nu=10_out=bs_pi=[PI]_samplesize=[SAMPLESIZE]_simnum=2_type=[SIMTYPE].csv`
`R=9_edge_mu=0.4_mu=[MU]_n_microbes=[K]_nu=10_out=main_effects_pi=[PI]_samplesize=[SAMPLESIZE]_simnum=2_type=[SIMTYPE].csv`
`R=9_edge_mu=0.4_mu=[MU]_n_microbes=[K]_nu=10_out=ms_pi=[PI]_samplesize=[SAMPLESIZE]_simnum=2_type=[SIMTYPE].csv`
`R=9_edge_mu=0.4_mu=[MU]_n_microbes=[K]_nu=10_out=xis_pi=[PI]_samplesize=[SAMPLESIZE]_simnum=2_type=[SIMTYPE].csv`
`R=9_edge_mu=0.4_mu=[MU]_n_microbes=[K]_nu=10_out=XYs_pi=[PI]_samplesize=[SAMPLESIZE]_simnum=2_type=[SIMTYPE].csv`

All simulations use random seed `734` for generating biologically relvant information, and the gaussian noise added to each sample is generated with different seeds for each simulation.

Gaussian seeds are generated with the following Julia code:

```
using Random

rng = MersenneTwister(734)
sample(rng, 100:9999, )
```

### Additional materials
The `additional_materials` directory contains information used for testing, as well as old implementations of [BayesianNetworkRegression.jl](https://github.com/samozm/BayesianNetworkRegression.jl).

`simulation/run_simulation.jl` is a julia script which will run the package.