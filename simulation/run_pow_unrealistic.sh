echo "unrealistic (power) sims"
#POWER
echo "mean = 0.8, pi = 0.0, r=$1, nu=$2, 8 microbes per sample, $3 samples"
julia simulation/run_simulation.jl -n 1 -m 0.8 -p 0.0 -r $1 -k 8 -s 2358 -u $2 -z $3

echo "mean = 0.8, pi = 0.0, r=$1, nu=$2, 15 microbes per sample, $3 samples"
julia simulation/run_simulation.jl -n 1 -m 0.8 -p 0.0 -r $1 -k 15 -s 2358 -u $2 -z $3

echo "mean = 0.8, pi = 0.0, r=$1, nu=$2, 22 microbes per sample, $3 samples"
julia simulation/run_simulation.jl -n 1 -m 0.8 -p 0.0 -r $1 -k 22 -s 2358 -u $2 -z $3

sleep 4s

echo "mean = 1.6, pi = 0.0, r=$1, nu=$2, 8 microbes per sample, $3 samples"
julia simulation/run_simulation.jl -n 1 -m 1.6 -p 0.0 -r $1 -k 8 -s 2358 -u $2 -z $3

echo "mean = 1.6, pi = 0.0, r=$1, nu=$2, 15 microbes per sample, $3 samples"
julia simulation/run_simulation.jl -n 1 -m 1.6 -p 0.0 -r $1 -k 15 -s 2358 -u $2 -z $3

echo "mean = 1.6, pi = 0.0, r=$1, nu=$2, 22 microbes per sample, $3 samples"
julia simulation/run_simulation.jl -n 1 -m 1.6 -p 0.0 -r $1 -k 22 -s 2358 -u $2 -z $3