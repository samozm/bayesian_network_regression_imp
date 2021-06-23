This folder contains scripts and output for JuliaCon 2021.
The simulations in the `data/` folder are obtained similarly to those described [here](/data/simulation/ReadMe.md), but only 30 total microbes are used, with t=8 and t=15 in each sample. Also, sparsity parameter (for the coefficient matrix) πₛ=0.3 is not used, only πₛ=0.1 and πₛ=0.8 are used.

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
bar(mses,series_annotations=round.(mses,digits=3),xaxis=(1:8),xlabel="case",ylabel="MSE",legend=false)
```
