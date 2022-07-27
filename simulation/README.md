This folder contains scripts to create the simulated data used for testing.

```
bash generate_unrealistic_sims.sh [SAMPLESIZE]
```
(to be run from the main folder `bayesian_network_regression_imp/`) generates simulation 1, which are "unrealistic". Simulated data is output in files in the `data/simulation/` folder. Currently 734 is used as the random seed for all simulations, but that can easily be changed by updating the file. Uses `unrealistic_sim.jl` to actually generate the data. \[SAMPLESIZE\] is a required parameter, the number of samples to generate.

```
bash generate_realistic_sims.sh [SAMPLESIZE]
```
(to be run from the main folder `bayesian_network_regression_imp/`) generates simulation 2, which are "realistic". This script generates simulated data for all 6 types of realistic simulations (additive phylogenetic, additive random, interaction phylogenetic, interaction random, functional redundancy phylogenetic, functional redundancy random - more information on these types of simulations is in notes.md) Simulated data is output in files in the `data/simulation/` folder. Currently 734 is used as the random seed for all simulations, but that can easily be changed by updating the file. Calls `generate_realistic_onetype` which in turn calls `realistic_sim.jl` to actually generate the data. \[SAMPLESIZE\] is a required parameter, the number of samples to generate. All edges are drawn from a normal distribution with mean 0.4 and standard deviation 1.

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


# Simulation

## Files
Each simulation consists of 4 files in the `data/simulation/[SIMTYPE]/` directory. Each file is named using a _-separated list of values used for the simulation (pi, mu, etc):
- The base file (out=XYs.csv) contains the X matrices and corresponding y vectors for the simulation.
- The b file (out=bs.csv) contains the true B matrix (column matrix) associated with the simulated data.
- The m file (out=ms.csv) contains the m value (which shows which taxa were "sampled" for each "sample" (row in the X matrix/entry in the y vector)). A 1 indicates that taxa was included in that sample, a 0 indicates that it wasn't.
- The xi file (out=xi) contais a vector of boolean values indicating whether each node is "turned on" - an edge has nonzero value if it connects two "turned on" nodes.

Realistic simulations also includes the following file:
- The main effects file (out=main_effects) contains the main effect for each node used to generate the results (since realistic simulations are generated using node main effects)

## Results
Results are in `results/simulation/[SIMTYPE]`.
For the model run, we will use 2358 as the random seed.

## Running Realistic simulations
We used `run_simulation_realistic.jl` to run realistic simulations, as follows:

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