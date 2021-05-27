# Simulation

### Files
Each simulation consists of 3 files:
- The base file (simulationX.csv) contains the X matrices and corresponding y vectors for the simulation.
- The b file (simulationX_bs.csv) contains the true B matrix (column matrix) associated with the simulated data.
- The m file (simulationX_ms.csv) contains the m value (which shows which taxa were "sampled" for each "sample" (row in the X matrix/entry in the y vector)). A 1 indicates that taxa was included in that sample, a 0 indicates that it wasn't.

In simulations 1 and 2, all taxa contribute to y (there are none with 0 contribution). Sample 1 has no interaction effect (B is a 0 column-matrix) and Sample 2 has a positive interaction effect (all entries in B are positive).


Questions:
1. Should the B-matrices for interaction actually be strictly positive/negative, or should they just be positive/negative definite?
