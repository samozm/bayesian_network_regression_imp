
echo "unrealistic simulations (sim 1)"

echo "case 1: mean = 0.8, pi = 0.1, samptaxa=20"
julia simulation/unrealistic_sim.jl -s 3323

echo "case 2: mean = 0.8, pi = 0.3, samptaxa=20"
julia simulation/unrealistic_sim.jl -s 3323 -p 0.3

echo "case 3: mean = 0.8, pi = 0.8, samptaxa=20"
julia simulation/unrealistic_sim.jl -s 3323 -p 0.8

echo "case 4: mean = 1.6, pi = 0.1, samptaxa=20"
julia simulation/unrealistic_sim.jl -s 3323 -m 1.6

echo "case 5: mean = 1.6, pi = 0.3, samptaxa=20"
julia simulation/unrealistic_sim.jl -s 3323 -m 1.6 -p 0.3

echo "case 6: mean = 1.6, pi = 0.8, samptaxa=20"
julia simulation/unrealistic_sim.jl -s 3323 -m 1.6 -p 0.8

echo "case 7: mean = 0.8, pi = 0.1, samptaxa=50"
julia simulation/unrealistic_sim.jl -s 3323 -k 50

echo "case 8: mean = 0.8, pi = 0.3, samptaxa=50"
julia simulation/unrealistic_sim.jl -s 3323 -p 0.3 -k 50

echo "case 9: mean = 0.8, pi = 0.8, samptaxa=50"
julia simulation/unrealistic_sim.jl -s 3323 -p 0.8 -k 50

echo "case 10: mean = 1.6, pi = 0.1, samptaxa=50"
julia simulation/unrealistic_sim.jl -s 3323 -m 1.6 -k 50

echo "case 11: mean = 1.6, pi = 0.3, samptaxa=50"
julia simulation/unrealistic_sim.jl -s 3323 -m 1.6 -p 0.3 -k 50

echo "case 12: mean = 1.6, pi = 0.8, samptaxa=50"
julia simulation/unrealistic_sim.jl -s 3323 -m 1.6 -p 0.8 -k 50
