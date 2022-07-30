This folder contains scripts to create the simulated data in Ozminkowski and Solis-Lemus (2022).

## Main scripts

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


## Complementary scripts

```
bash generate_realistic_onetype.sh [TYPE] [SAMPLES]
```
(to be run from the main folder `bayesian_network_regression_imp/`) generates one type of realistic simulation (of the 6 types described above). It runs the model on all combinations of pi=0.3,0.8 (density parameter of the coefficient matrix), mu=0.8,1.6 (magnitude of node coefficients), and microbes per sample=8,22 (density of adjacency matrix). 

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
bash run_sim1.sh
```
(to be run from the main folder `bayesian_network_regression_imp/`) calls `simulation/run_unrealistic.sh` four times, with R=5 nu=10, R=9 nu=10, R=10 nu=15, and R=15 nu=20, plus `simulation/run_pow_unrealistic.sh` with R=9 nu=10 for power simulations (B=**0**). 100 samples is used for all runs.

```
bash run_pow_unrealistic.sh [R] [NU] [SAMPLESIZE]
```
(to be run from the main folder `bayesian_network_regression_imp/`) takes as input R value, nu value, and number of samples where R is the dimension of the latent variable u used in the model. It runs the model on all combinations of mu=0.8,1.6 (magnitude of edge coefficients), and microbes per sample=8,15,22 (density of adjacency matrix) with pi=0 (resulting in B=**0**). 

```
bash run_Xsamp_sim1.sh [SAMPLESIZE]
```
(to be run from the main folder `bayesian_network_regression_imp/`) calls `simulation/run_unrealistic` with R=9 nu=10 and sample size \[SAMPLESIZE\].

```
bash run_sim2.sh
```
(to be run from the main folder `bayesian_network_regression_imp/`) calls `simulation/run_realistic` six times, with R=9 nu=10, once for each type of realistic simulation (listed above, more information in notes.md) 100 samples is used for all runs.

```
bash run_Xsamp_sim2.sh [SAMPLESIZE]
```
(to be run from the main folder `bayesian_network_regression_imp/`) calls `simulation/run_realistic` with R=9 nu=10 and sample size \[SAMPLESIZE\] for each type of realistic simulation.

```
julia run_simulation.jl
```
(to be run from the main folder `bayesian_network_regression_imp/`) handles actually running the Bayesian Network Regression and outputting the results (in `results/simulation/[SIMTYPE]`, where \[SIMTYPE\] is either "realistic" or "unrealistic") into four files: one contains MSE values, one contains posterior probabilities of influence for each node (as well as whether they were truly influential), one contains the following posterior values for each edge coefficient: 2.5th %ile, 5th %ile, mean, 95%ile, and 97.5%ile (for forming credible intervals) and one contains posterior values for mu, the overall mean. Each output file also contains the values of pi, mu, nu, R, and number of microbes used for that run. See file for argument options.