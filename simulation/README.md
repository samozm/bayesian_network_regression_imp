This folder contains scripts to create the simulated data used for testing.

```
bash generate_unrealistic_sims.sh
```
(to be run from the main folder `bayesian_network_regression_imp/`) generates simulation 1 (cases 1-24), which are "unrealistic". Simulated data is output in files in the `data/simulation/` folder. Currently 3323 is used as the random seed for all simulations, but that can easily be changed by updating the file. 

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

    1. using ape, simulate 1 phylogenetic tree with 100 total microbes
    2. calculate the phylogenetic distance between each pair of microbes (again using ape)
    3. for each sample (70 total) randomly (uniform) select t microbes
        - t = 20, 50
        - for A[k,l] if k and l are both chosen A[k,l] is the inverse of the phylogenetic distance between microbes k and l. Otherwise A[k,l] = 0


To generate **B**<sub>0</sub>:

    1. Generate 100 binary indicators ξ₁⁰,..., ξₜ⁰ independently from Ber(πₛ) (indicator of influential node)
        - πₛ = 0.1, 0.3, 0.8
    2. If ξₖ⁰ = ξₗ⁰ = 1$ generate the edge coefficient between k and l from N(μₛ,1)
        - μₛ = 0.8,1.6
       Else set the edge coefficient between k and l to 0


Questions:
1. For realistic: should the B-matrices for interaction actually be strictly positive/negative, or should they just be positive/negative definite?

