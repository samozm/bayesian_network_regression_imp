
echo "juliacon unrealistic simulations (sim 1)"

echo "case 1: mean = 0.8, pi = 0.1, samptaxa=10"
julia simulation/unrealistic_sim.jl -s 3323 -c 1 -k 10 -t 40 -j

echo "case 2: mean = 0.8, pi = 0.8, samptaxa=10"
julia simulation/unrealistic_sim.jl -s 3323 -c 2 -p 0.8 -k 10 -t 40 -j

echo "case 3: mean = 1.6, pi = 0.1, samptaxa=10"
julia simulation/unrealistic_sim.jl -s 3323 -c 3 -m 1.6 -k 10 -t 40 -j

echo "case 4: mean = 1.6, pi = 0.8, samptaxa=10"
julia simulation/unrealistic_sim.jl -s 3323 -c 4 -m 1.6 -p 0.8 -k 10 -t 40 -j

echo "case 5: mean = 0.8, pi = 0.1, samptaxa=20"
julia simulation/unrealistic_sim.jl -s 3323 -c 5 -k 20 -t 40 -j

echo "case 6: mean = 0.8, pi = 0.8, samptaxa=50"
julia simulation/unrealistic_sim.jl -s 3323 -c 6 -p 0.8 -k 20 -t 40 -j

echo "case 7: mean = 1.6, pi = 0.1, samptaxa=50"
julia simulation/unrealistic_sim.jl -s 3323 -c 7 -m 1.6 -k 20 -t 40 -j

echo "case 8: mean = 1.6, pi = 0.8, samptaxa=50"
julia simulation/unrealistic_sim.jl -s 3323 -c 8 -m 1.6 -p 0.8 -k 20 -t 40 -j
