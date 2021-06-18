# Simulation

## Files
Each simulation consists of 3 files:
- The base file (simulationX.csv) contains the X matrices and corresponding y vectors for the simulation.
- The b file (simulationX_bs.csv) contains the true B matrix (column matrix) associated with the simulated data.
- The m file (simulationX_ms.csv) contains the m value (which shows which taxa were "sampled" for each "sample" (row in the X matrix/entry in the y vector)). A 1 indicates that taxa was included in that sample, a 0 indicates that it wasn't.

Simulation 1 (cases 1-6) is "unrealistic". Generated as follows:

y<sub>i</sub> &= $\mu$<sub>0</sub> + <**A**<sub>i</sub>,**B**<sub>0</sub>><sub>F</sub> + $\epsilon$<sub>i</sub>

$\epsilon$<sub>i</sub> ~ N(0,$\tau$<sub>0</sub><sup>2</sup>)

$\tau$<sub>0</sub><sup>2</sup> &= 1


To generate $\textbf{A}_i$:

    1. using ape, simulate 70 phylogenetic trees with
        - 20 (?) total microbes
        - 10 (?) microbes per tree
    2. calculate the phylogenetic distance between each pair of microbes (again using ape)


To generate **B**<sub>0</sub>:

    1. Generate 20 binary indicators $\xi$<sub>1</sub><sup>0</sup>,\hdots, $\xi$ _{20}^0$ independently from Ber($\pi_{s}$) (indicator of influential node)
        - $\pi_s$ = 0.1, 0.3, 0.8
    2. If $\xi_k^0 = \xi_l^0 = 1$ generate the edge coefficient between k and l from $N(\mu_s,1)$
        - $\mu_s$ = 0.8,1.6
       Else set the edge coefficient between k and l to 0


Questions:
1. For realistic: should the B-matrices for interaction actually be strictly positive/negative, or should they just be positive/negative definite?
