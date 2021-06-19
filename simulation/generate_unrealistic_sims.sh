
echo "unrealistic simulations (sim 1)"

echo "case 1: mean = 0.8, pi = 0.1"
julia simulation/unrealistic_sim.jl -s 3323 -c 1

echo "case 2: mean = 0.8, pi = 0.3"
julia simulation/unrealistic_sim.jl -s 3323 -c 2 -p 0.3

echo "case 3: mean = 0.8, pi = 0.8"
julia simulation/unrealistic_sim.jl -s 3323 -c 3 -p 0.8

echo "case 4: mean = 1.6, pi = 0.1"
julia simulation/unrealistic_sim.jl -s 3323 -c 4 -m 1.6

echo "case 5: mean = 1.6, pi = 0.3"
julia simulation/unrealistic_sim.jl -s 3323 -c 5 -m 1.6 -p 0.3

echo "case 6: mean = 1.6, pi = 0.8"
julia simulation/unrealistic_sim.jl -s 3323 -c 6 -m 1.6 -p 0.8
