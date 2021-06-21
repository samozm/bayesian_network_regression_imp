
echo "unrealistic simulations (sim 1)"

echo "case 1: mean = 0.8, pi = 0.1"
julia simulation/run_simulation.jl -n 1 -c 1

echo "case 2: mean = 0.8, pi = 0.3"
julia simulation/run_simulation.jl -n 1 -c 2 

echo "case 3: mean = 0.8, pi = 0.8"
julia simulation/run_simulation.jl -n 1 -c 3

echo "case 4: mean = 1.6, pi = 0.1"
julia simulation/run_simulation.jl -n 1 -c 4

echo "case 5: mean = 1.6, pi = 0.3"
julia simulation/run_simulation.jl -n 1 -c 5

echo "case 6: mean = 1.6, pi = 0.8"
julia simulation/run_simulation.jl -n 1 -c 6
