This folder contains scripts to create the simulated data used for testing.

```
bash generate_unrealistic_sims.sh
```
(to be run from the main folder `bayesian_network_regression_imp/`) generates simulations 1-6, which are "unrealistic". For an in-depth explanation of how they're generated, see the README for the data folder [here](/data/simulation/ReadMe.md). Simulated data is output in files in the `data/simulation/` folder. Currently 3323 is used as the random seed for all simulations, but that can easily be changed by updating the file. 
