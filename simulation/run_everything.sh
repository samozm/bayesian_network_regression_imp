echo "running EVERYTHING - this will take ~ 5-6 days"

julia --optimize=0 --math-mode=ieee --check-bounds=yes bayesian_network_regression_imp/simulation/run_simulation_unrealistic.jl

julia --optimize=0 --math-mode=ieee --check-bounds=yes bayesian_network_regression_imp/simulation/run_simulation_realistic.jl