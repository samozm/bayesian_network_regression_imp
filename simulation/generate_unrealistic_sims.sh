
echo "unrealistic simulations\n"

echo "sim 1: mean = 0.8, pi = 0.1"
julia simulation/unrealistic_sim.jl -s 3323 -i 1

echo "sim 2: mean = 0.8, pi = 0.3"
julia simulation/unrealistic_sim.jl -s 3323 -i 2 -p 0.3

echo "sim 3: mean = 0.8, pi = 0.8"
julia simulation/unrealistic_sim.jl -s 3323 -i 3 -p 0.8

echo "sim 4: mean = 1.6, pi = 0.1"
julia simulation/unrealistic_sim.jl -s 3323 -i 4 -m 1.6

echo "sim 5: mean = 1.6, pi = 0.3"
julia simulation/unrealistic_sim.jl -s 3323 -i 5 -m 1.6 -p 0.3

echo "sim 6: mean = 1.6, pi = 0.8"
julia simulation/unrealistic_sim.jl -s 3323 -i 6 -m 1.6 -p 0.8
