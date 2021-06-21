
echo "unrealistic simulations for juliacon (sim 1)"

echo "case 1: mean = 0.8, pi = 0.1, 10 microbes per sample"
julia simulation/run_simulation.jl -n 1 -c 1 -j

echo "case 2: mean = 0.8, pi = 0.8, 10 microbes per sample"
julia simulation/run_simulation.jl -n 1 -c 2 -j

echo "case 3: mean = 1.6, pi = 0.1, 10 microbes per sample"
julia simulation/run_simulation.jl -n 1 -c 3 -j

echo "case 4: mean = 1.6, pi = 0.8, 10 microbes per sample"
julia simulation/run_simulation.jl -n 1 -c 4 -j

echo "case 5: mean = 0.8, pi = 0.1, 10 microbes per sample"
julia simulation/run_simulation.jl -n 1 -c 5 -j

echo "case 6: mean = 0.8, pi = 0.8, 10 microbes per sample"
julia simulation/run_simulation.jl -n 1 -c 6 -j

echo "case 7: mean = 1.6, pi = 0.1, 10 microbes per sample"
julia simulation/run_simulation.jl -n 1 -c 7 -j

echo "case 8: mean = 1.6, pi = 0.8, 10 microbes per sample"
julia simulation/run_simulation.jl -n 1 -c 8 -j
