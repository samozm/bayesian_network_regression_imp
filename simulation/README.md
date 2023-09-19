This folder contains scripts to create the simulated data in Ozminkowski and Solis-Lemus (2022).

## Main scripts
```
bash generate_all_sims.jl
```
(to be run from the main folder `bayesian_network_regression_imp/`) generates all types of unrealistic and realistic simulations (of the 6 types described below). It runs the model on all combinations of pi=0.3,0.8 (density parameter of the coefficient matrix), mu=0.8,1.6 (magnitude of node coefficients), microbes per sample=8,22 (density of adjacency matrix), samplesize (=100,500 for theoretical, 500,1000 for realistic), and simulation type (additive phylogenetic, additive random, interaction phylogenetic, interaction random, functional redundancy phylogenetic, functional redundancy random -- for realistic only). 


## Complementary scripts


`sim_trees.R` contains source code used to simulate trees with `ape`.

Unrealistic (theoretical) simulations:
```
bash generate_unrealistic_sims.sh [SAMPLESIZE]
```

Note that this script should be run from the main folder `bayesian_network_regression_imp/`, and it generates simulation 1, which are "unrealistic". Simulated data is output in files in the `data/simulation/` folder. Currently 734 is used as the random seed for all simulations, but that can easily be changed by updating the file. The script `generate_unrealistic_sims.sh` calls `unrealistic_sim.jl` to actually generate the data. \[SAMPLESIZE\] is a required parameter, the number of samples to generate.

Realistic simulations:
```
bash generate_realistic_sims.sh [SAMPLESIZE]
```

Note that this script should be run from the main folder  `bayesian_network_regression_imp/`), and it generates simulation 2, which are "realistic". This script generates simulated data for all 6 types of realistic simulations (additive phylogenetic, additive random, interaction phylogenetic, interaction random, functional redundancy phylogenetic, functional redundancy random). Simulated data is output in files in the `data/simulation/` folder. Currently 734 is used as the random seed for all simulations, but that can easily be changed by updating the file. The script `generate_realistic_sims.sh` calls `generate_realistic_onetype` which in turn calls `realistic_sim.jl` to actually generate the data. \[SAMPLESIZE\] is a required parameter, the number of samples to generate. All edges are drawn from a normal distribution with mean 0.4 and standard deviation 1.



```
julia simulation/run_simulation_unrealistic.jl
```
(to be run from the main folder `bayesian_network_regression_imp/`). It runs the model on all combinations of *pi=0.3,0.8* (density parameter of the coefficient matrix), *mu=0.8,1.6* (magnitude of edge coefficients),  *microbes per sample=8,15,22* (density of adjacency matrix), *samplesize=100,500*, and latent dimension of u *R=5,7,9* for the unrealistic simulations. 

```
julia simulation/run_simulation_unrealistic.jl
```
(to be run from the main folder `bayesian_network_regression_imp/`). It runs the model on all combinations of *pi=0.3,0.8* (density parameter of the coefficient matrix), *mu=0.8,1.6* (magnitude of edge coefficients), *microbes per sample=8,22* (density of adjacency matrix), *samplesize=500,1000*, and *simtype=additive_phylo, additive_random, interaction_phylo, interaction_random, redundant_phylo, redundant_random* for the realistic simulations. 


## Running Realistic simulations
We use `run_simulation_realistic.jl` to run realistic simulations, as follows:

```
julia --optimize=0 --math-mode=ieee --check-bounds=yes simulation/run_simulation_realistic.jl
```

## Running Unrealistic simulations
For the unrealistic simulations `simulation/run_simulation_unrealistic.jl` is run
```
julia --optimize=0 --math-mode=ieee --check-bounds=yes simulation/run_simulation_unrealistic.jl
```


## Additional (old) scripts


```
julia run_simulation.jl
```
(to be run from the main folder `bayesian_network_regression_imp/`) handles actually running the Bayesian Network Regression and outputting the results (in `results/simulation/[SIMTYPE]`, where \[SIMTYPE\] is either "realistic" or "unrealistic") into four files: one contains MSE values, one contains posterior probabilities of influence for each node (as well as whether they were truly influential), one contains the following posterior values for each edge coefficient: 2.5th %ile, 5th %ile, mean, 95%ile, and 97.5%ile (for forming credible intervals) and one contains posterior values for mu, the overall mean. Each output file also contains the values of pi, mu, nu, R, and number of microbes used for that run. See file for argument options.