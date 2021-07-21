This folder contains scripts to create the simulated data used for testing.

```
bash generate_unrealistic_sims.sh
```
(to be run from the main folder `bayesian_network_regression_imp/`) generates simulation 1 (cases 1-24), which are "unrealistic". Simulated data is output in files in the `data/simulation/` folder. Currently 734 is used as the random seed for all simulations, but that can easily be changed by updating the file. Uses `unrealistic_sim.jl` to actually generate the data.

```
bash run_unrealistic.sh X Y 
```
(to be run from the main folder `bayesian_network_regression_imp/`) takes as input R value (X) and nu value (Y) where R is the dimension of the latent variable u used in the model. It runs the model on all combinations of pi=0.3,0.8 (density parameter of the coefficient matrix), mu=0.8,1.6 (magnitude of edge coefficients), and microbes per sample=8,15,22 (density of adjacency matrix). 

```
bash run_sim1.sh
```
(to be run from the main folder `bayesian_network_regression_imp/`) calls `simulation/run_unrealistic` three times, with R=5 nu=10, R=10 nu=15, and R=15 nu=20.

```
julia run_simulation.jl
```
(to be run from the main folder `bayesian_network_regression_imp/`) handles actually running the Bayesian Network Regression and outputting the results (in `results/simulation/`) into 3 files: one contains MSE values, one contains posterior probabilities of influence for each node (as well as whether they were truly influential), and one contains the following posterior values for each edge coefficient: 2.5th %ile, 5th %ile, mean, 95%ile, and 97.5%ile (for forming credible intervals). Each output file also contains the values of pi, mu, nu, R, and number of microbes used for that run. See file for argument options.

# Simulation

## Files
Each simulation consists of 3 files:
- The base file (simulationX.csv) contains the X matrices and corresponding y vectors for the simulation.
- The b file (simulationX_bs.csv) contains the true B matrix (column matrix) associated with the simulated data.
- The m file (simulationX_ms.csv) contains the m value (which shows which taxa were "sampled" for each "sample" (row in the X matrix/entry in the y vector)). A 1 indicates that taxa was included in that sample, a 0 indicates that it wasn't.

Simulation 1 (cases 1-6) is "unrealistic". Generated as follows:

y<sub>i</sub> = μ<sub>0</sub> + <**A**<sub>i</sub>,**B**<sub>0</sub>><sub>F</sub> + ϵ <sub>i</sub>

ϵ<sub>i</sub> ~ N(0,τ<sub>0</sub><sup>2</sup>)

τ<sub>0</sub><sup>2</sup> = 1


To generate **A**<sub>i</sub>:

    1. using ape, simulate 1 phylogenetic tree with 30 total microbes
    2. calculate the phylogenetic distance between each pair of microbes (again using ape)
    3. for each sample (100 total) randomly (uniform) select t microbes
        - t = 8,15,22
        - for A[k,l] if k and l are both chosen A[k,l] is the inverse of the phylogenetic distance between microbes k and l. Otherwise A[k,l] = 0


To generate **B**<sub>0</sub>:

    1. Generate 100 binary indicators ξ₁⁰,..., ξₜ⁰ independently from Ber(πₛ) (indicator of influential node)
        - πₛ = 0.3, 0.8
    2. If ξₖ⁰ = ξₗ⁰ = 1$ generate the edge coefficient between k and l from N(μₛ,1)
        - μₛ = 0.8,1.6
       Else set the edge coefficient between k and l to 0

## Results
We will run each of the simulated datasets through the model three times, once using R=5 and nu=10, once using R=10 and nu=15, and once using R=15 and nu=20.
For the model run, we will use 2358 as the random seed.