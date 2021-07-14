This folder contains scripts and output for JuliaCon 2021.

ALL SCRIPTS SHOULD BE RUN FROM THE BASE FOLDER (`bayesian_network_regression_imp/`). FAILURE TO DO SO COULD CAUSE UNEXPECTED BEHAVIOR.

The simulations in the `data/` folder are obtained similarly to those described [here](/data/simulation/ReadMe.md), but only 30 total microbes are used, with t=8 and t=15 in each sample. Also, sparsity parameter (for the coefficient matrix) πₛ=0.3 is not used, only πₛ=0.1 and πₛ=0.8 are used. For generating the data, seed 734 was used.

Results are obtained by first running
```
zsh juliacon/src/generate_unrealistic_data.sh
```
to generate the data. Then we run
```
zsh juliacon/src/run_sim1.sh
```
to execute the model. Results (true ξ values, posterior means for B/γ).

To plot MSEs add them to a vector called `mse` and run in julia:
```
bar(mse,series_annotations=round.(mse,digits=3),xaxis=(1:8),xlabel="case",ylabel="MSE",legend=false)
```

MSEs for each case are:

 Case |   1   |   2   |   3   |   4   |   5   |   6   |   7   |   8   
------|-------|-------|-------|-------|-------|-------|-------|-------
 MSE  | 0.046 | 0.553 | 0.014 | 0.849 | 0.011 | 0.441 | 0.006 | 1.421


To produce plots for the gammas run in julia:
```
include("juliacon/src/analysis.jl")

# parameters are simulation number and case number
analyze_γ(1,1)
analyze_γ(1,2)
analyze_γ(1,3)
analyze_γ(1,4)
analyze_γ(1,5)
analyze_γ(1,6)
analyze_γ(1,7)
analyze_γ(1,8)
```
