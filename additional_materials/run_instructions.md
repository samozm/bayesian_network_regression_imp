
1. Download BayesianNetworkRegression from here: [https://github.com/samozm/BayesianNetworkRegression.jl](https://github.com/samozm/BayesianNetworkRegression.jl)

2. Grab the files `run_simulation.jl`, `prepare_results.jl`, and `realistic/data/edge_mu=0.4_mu=1.6_n_microbes=8_out=[OUT]_pi=0.3_samplesize=500_simnum=2_type=redundant_phylo.csv` from `submit2.chtc.wisc.edu/home/ozminkowski` into your current directory, where OUT is each of the following (XYs, bs, xis) (but preserve the directory structure for the data files, so `current_directory/realistic/data/...`)

3. Also move `run_and_prep.sh` into the current directory

4. create a new directory called `realistic-results` to hold the results

5. create a new directory to hold julia enviroment (ex. `my_project/`)

6. start julia and enter pkg mode
```
julia --project=my_project 
]
```

7. add the following julia packages
```
add LinearAlgebra,CSV,ArgParse,Random, DataFrames, StatsBase, InvertedIndices, ProgressMeter, Distributions, StaticArrays,TypedTables,DrWatson,MCMCDiagnosticTools,JLD2
```

8. add the development version of BayesianNetworkRegression
```
dev BayesianNetworkRegression.jl
```
(with the path to BayesianNetworkRegression.jl if it's not in the current directory)

9. exit julia and run the following from the command line

```
bash run_and_prep.sh 1.6 0.3 8 500 redundant_phylo
```

10. results will be in files called `realistic-results/R=9_edge_mu=0.4_mu=1.6_n_microbes=8_nu=10_out=[OUT]_pi=0.3_samplesize=500_simnum=2_type=redundant_phylo.csv` where out is each of (time,nodes,edges,mu,MSE,psrf)