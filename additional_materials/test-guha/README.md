### Testing
Tests are run on Simulation 1 Case 1, Simulation 2 cases 1 and 2 as described in [Guha and Rodriguez (2018)](https://www.tandfonline.com/doi/full/10.1080/01621459.2020.1772079).

To run all 3 cases in order run:
```
bash additional_materials/test/run_all.sh --nburn X --nsamp Y
```

Where X and Y are integers, X is the number of burn-in gibbs samples to discard
and Y is the number of gibbs samples to use. For testing, we used 30000 burn-in samples
and 20000 post-burn samples (the same number used by Guha &amp; Rodriguez).

This runs my implementation as well as Guha's on the same data.

Final ξ values (posterior probability of a node being selected) and MSE for
my implementation and Guha's implementation are output to text.

Graphs of γ values (coefficients for each edge) are also created in `plots/test/` (Guha's in red, ours in blue) and automatically opened. 
These tests appear to be consistent between our implementation and Guha's.