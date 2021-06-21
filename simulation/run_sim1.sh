
echo "unrealistic simulations (sim 1)"

echo "case 1: mean = 0.8, pi = 0.1, 20 microbes per sample"
julia simulation/run_simulation.jl -n 1 -c 1

echo "case 2: mean = 0.8, pi = 0.3, 20 microbes per sample"
julia simulation/run_simulation.jl -n 1 -c 2

echo "case 3: mean = 0.8, pi = 0.8, 20 microbes per sample"
julia simulation/run_simulation.jl -n 1 -c 3

echo "case 4: mean = 1.6, pi = 0.1, 20 microbes per sample"
julia simulation/run_simulation.jl -n 1 -c 4

echo "case 5: mean = 1.6, pi = 0.3, 20 microbes per sample"
julia simulation/run_simulation.jl -n 1 -c 5

echo "case 6: mean = 1.6, pi = 0.8, 20 microbes per sample"
julia simulation/run_simulation.jl -n 1 -c 6

echo "case 7: mean = 0.8, pi = 0.1, 50 microbes per sample"
julia simulation/run_simulation.jl -n 1 -c 7

echo "case 8: mean = 0.8, pi = 0.3, 50 microbes per sample"
julia simulation/run_simulation.jl -n 1 -c 8

echo "case 9: mean = 0.8, pi = 0.8, 50 microbes per sample"
julia simulation/run_simulation.jl -n 1 -c 9

echo "case 10: mean = 1.6, pi = 0.1, 50 microbes per sample"
julia simulation/run_simulation.jl -n 1 -c 10

echo "case 11: mean = 1.6, pi = 0.3, 50 microbes per sample"
julia simulation/run_simulation.jl -n 1 -c 11

echo "case 12: mean = 1.6, pi = 0.8, 50 microbes per sample"
julia simulation/run_simulation.jl -n 1 -c 12
