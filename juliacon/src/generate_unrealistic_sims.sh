
echo "juliacon unrealistic simulations (sim 1) - 30 nodes total"

echo "case 1: mean = 0.8, pi = 0.1, samptaxa=8"
julia simulation/unrealistic_sim.jl -s 3323 -c 1 -k 8 -t 30 -j

echo "case 2: mean = 0.8, pi = 0.8, samptaxa=8"
julia simulation/unrealistic_sim.jl -s 3323 -c 2 -p 0.8 -k 8 -t 30 -j

echo "case 3: mean = 1.6, pi = 0.1, samptaxa=8"
julia simulation/unrealistic_sim.jl -s 3323 -c 3 -m 1.6 -k 8 -t 30 -j

echo "case 4: mean = 1.6, pi = 0.8, samptaxa=8"
julia simulation/unrealistic_sim.jl -s 3323 -c 4 -m 1.6 -p 0.8 -k 8 -t 30 -j

echo "case 5: mean = 0.8, pi = 0.1, samptaxa=15"
julia simulation/unrealistic_sim.jl -s 3323 -c 5 -k 15 -t 30 -j

echo "case 6: mean = 0.8, pi = 0.8, samptaxa=15"
julia simulation/unrealistic_sim.jl -s 3323 -c 6 -p 0.8 -k 15 -t 30 -j

echo "case 7: mean = 1.6, pi = 0.1, samptaxa=15"
julia simulation/unrealistic_sim.jl -s 3323 -c 7 -m 1.6 -k 15 -t 30 -j

echo "case 8: mean = 1.6, pi = 0.8, samptaxa=15"
julia simulation/unrealistic_sim.jl -s 3323 -c 8 -m 1.6 -p 0.8 -k 15 -t 30 -j
