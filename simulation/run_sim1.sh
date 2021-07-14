
echo "unrealistic simulations (sim 1)"

echo "case 1: mean = 0.8, pi = 0.1, r=5, 20 microbes per sample"
julia simulation/run_simulation.jl -n 1 -m 0.8 -p 0.1 -r 5 -k 20 -s 734 

echo "case 2: mean = 0.8, pi = 0.3, r=5, 20 microbes per sample"
julia simulation/run_simulation.jl -n 1 -m 0.8 -p 0.3 -r 5 -k 20 -s 734 

echo "case 3: mean = 0.8, pi = 0.8, r=5, 20 microbes per sample"
julia simulation/run_simulation.jl -n 1 -m 0.8 -p 0.8 -r 5 -k 20 -s 734 

echo "case 4: mean = 1.6, pi = 0.1, r=5, 20 microbes per sample"
julia simulation/run_simulation.jl -n 1 -m 1.6 -p 0.1 -r 5 -k 20 -s 734 

echo "case 5: mean = 1.6, pi = 0.3, r=5, 20 microbes per sample"
julia simulation/run_simulation.jl -n 1 -m 1.6 -p 0.3 -r 5 -k 20 -s 734 

echo "case 6: mean = 1.6, pi = 0.8, r=5, 20 microbes per sample"
julia simulation/run_simulation.jl -n 1 -m 1.6 -p 0.8 -r 5 -k 20 -s 734 

echo "case 7: mean = 0.8, pi = 0.1, r=5, 50 microbes per sample"
julia simulation/run_simulation.jl -n 1 -m 0.8 -p 0.1 -r 5 -k 50 -s 734 

echo "case 8: mean = 0.8, pi = 0.3, r=5, 50 microbes per sample"
julia simulation/run_simulation.jl -n 1 -m 0.8 -p 0.3 -r 5 -k 50 -s 734 

echo "case 9: mean = 0.8, pi = 0.8, r=5, 50 microbes per sample"
julia simulation/run_simulation.jl -n 1 -m 0.8 -p 0.8 -r 5 -k 50 -s 734 

echo "case 10: mean = 1.6, pi = 0.1, r=5, 50 microbes per sample"
julia simulation/run_simulation.jl -n 1 -m 1.6 -p 0.1 -r 5 -k 50 -s 734 

echo "case 11: mean = 1.6, pi = 0.3, r=5, 50 microbes per sample"
julia simulation/run_simulation.jl -n 1 -m 1.6 -p 0.3 -r 5 -k 50 -s 734 

echo "case 12: mean = 1.6, pi = 0.8, r=5, 50 microbes per sample"
julia simulation/run_simulation.jl -n 1 -m 1.6 -p 0.8 -r 5 -k 50 -s 734 


