
echo "unrealistic simulations for juliacon (sim 1)"

echo "case 1: mean = 0.8, pi = 0.1, 8 microbes per sample"
julia simulation/run_simulation.jl -n 1 -m 0.8 -p 0.1 -r 5 -k 8 -j -s 734 

echo "case 2: mean = 0.8, pi = 0.8, 8 microbes per sample"
julia simulation/run_simulation.jl -n 1 -m 0.8 -p 0.8 -r 5 -k 8 -j -s 734 

echo "case 3: mean = 1.6, pi = 0.1, 8 microbes per sample"
julia simulation/run_simulation.jl -n 1 -m 1.6 -p 0.1 -r 5 -k 8 -j -s 734 

echo "case 4: mean = 1.6, pi = 0.8, 8 microbes per sample"
julia simulation/run_simulation.jl -n 1 -m 1.6 -p 0.8 -r 5 -k 8 -j -s 734 

echo "case 5: mean = 0.8, pi = 0.1, 15 microbes per sample"
julia simulation/run_simulation.jl -n 1  -m 0.8 -p 0.1 -r 5 -k 15 -j -s 734 

echo "case 6: mean = 0.8, pi = 0.8, 15 microbes per sample"
julia simulation/run_simulation.jl -n 1  -m 0.8 -p 0.8 -r 5 -k 15 -j -s 734 

echo "case 7: mean = 1.6, pi = 0.1, 15 microbes per sample"
julia simulation/run_simulation.jl -n 1  -m 1.6 -p 0.1 -r 5 -k 15 -j -s 734 

echo "case 8: mean = 1.6, pi = 0.8, 15 microbes per sample"
julia simulation/run_simulation.jl -n 1  -m 1.6 -p 0.8 -r 5 -k 15 -j -s 734 
