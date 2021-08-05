
echo "unrealistic simulations (sim 1)"

echo "mean = 0.8, pi = 0.3, samptaxa=8, $1 samples"
julia simulation/unrealistic_sim.jl -s 734 -m 0.8 -p 0.3 -t 30 -k 8 --samplesize $1

echo "mean = 0.8, pi = 0.3, samptaxa=15, $1 samples"
julia simulation/unrealistic_sim.jl -s 734 -m 0.8 -p 0.3 -t 30 -k 15 --samplesize $1

echo "mean = 0.8, pi = 0.3, samptaxa=22, $1 samples"
julia simulation/unrealistic_sim.jl -s 734 -m 0.8 -p 0.3 -t 30 -k 22 --samplesize $1

echo "mean = 0.8, pi = 0.8, samptaxa=8, $1 samples"
julia simulation/unrealistic_sim.jl -s 734 -m 0.8 -p 0.8 -t 30 -k 8 --samplesize $1

echo "mean = 0.8, pi = 0.8, samptaxa=15, $1 samples"
julia simulation/unrealistic_sim.jl -s 734 -m 0.8 -p 0.8 -t 30 -k 15 --samplesize $1

echo "mean = 0.8, pi = 0.8, samptaxa=22, $1 samples"
julia simulation/unrealistic_sim.jl -s 734 -m 0.8 -p 0.8 -t 30 -k 22 --samplesize $1

echo "mean = 1.6, pi = 0.3, samptaxa=8, $1 samples"
julia simulation/unrealistic_sim.jl -s 734 -m 1.6 -p 0.3 -t 30 -k 8 --samplesize $1

echo "mean = 1.6, pi = 0.3, samptaxa=15, $1 samples"
julia simulation/unrealistic_sim.jl -s 734 -m 1.6 -p 0.3 -t 30 -k 15 --samplesize $1

echo "mean = 1.6, pi = 0.3, samptaxa=22, $1 samples"
julia simulation/unrealistic_sim.jl -s 734 -m 1.6 -p 0.3 -t 30 -k 22 --samplesize $1

echo "mean = 1.6, pi = 0.8, samptaxa=8, $1 samples"
julia simulation/unrealistic_sim.jl -s 734 -m 1.6 -p 0.8 -t 30 -k 8 --samplesize $1

echo "mean = 1.6, pi = 0.8, samptaxa=15, $1 samples"
julia simulation/unrealistic_sim.jl -s 734 -m 1.6 -p 0.8 -t 30 -k 15 --samplesize $1

echo "mean = 1.6, pi = 0.8, samptaxa=22, $1 samples"
julia simulation/unrealistic_sim.jl -s 734 -m 1.6 -p 0.8 -t 30 -k 22 --samplesize $1


# POWER
echo "POWER"

echo "mean = 0.8, pi = 0.0, samptaxa=8, $1 samples"
julia simulation/unrealistic_sim.jl -s 734 -m 0.8 -p 0.0 -t 30 -k 8 --samplesize $1

echo "mean = 0.8, pi = 0.0, samptaxa=15, $1 samples"
julia simulation/unrealistic_sim.jl -s 734 -m 0.8 -p 0.0 -t 30 -k 15 --samplesize $1

echo "mean = 0.8, pi = 0.0, samptaxa=22, $1 samples"
julia simulation/unrealistic_sim.jl -s 734 -m 0.8 -p 0.0 -t 30 -k 22 --samplesize $1

echo "mean = 1.6, pi = 0.0, samptaxa=8, $1 samples"
julia simulation/unrealistic_sim.jl -s 734 -m 1.6 -p 0.0 -t 30 -k 8 --samplesize $1

echo "mean = 1.6, pi = 0.0, samptaxa=15, $1 samples"
julia simulation/unrealistic_sim.jl -s 734 -m 1.6 -p 0.0 -t 30 -k 15 --samplesize $1

echo "mean = 1.6, pi = 0.0, samptaxa=22, $1 samples"
julia simulation/unrealistic_sim.jl -s 734 -m 1.6 -p 0.0 -t 30 -k 22 --samplesize $1
