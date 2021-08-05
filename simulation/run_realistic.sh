#!/bin/bash

echo "mean = 0.8, pi = 0.3, r=$1, nu=$2, $4 samples, 8 microbes per sample"
julia simulation/run_simulation.jl -n 2 -m 0.8 -p 0.3 -r $1 -k 8 -s 2358 -u $2 --simtype $3 -e 0.4 -z $4

echo "mean = 0.8, pi = 0.3, r=$1, nu=$2, $4 samples, 22 microbes per sample"
julia simulation/run_simulation.jl -n 2 -m 0.8 -p 0.3 -r $1 -k 22 -s 2358 -u $2 --simtype $3 -e 0.4 -z $4


sleep 4s

echo "mean = 0.8, pi = 0.8, r=$1, nu=$2, $4 samples, 8 microbes per sample"
julia simulation/run_simulation.jl -n 2 -m 0.8 -p 0.8 -r $1 -k 8 -s 2358 -u $2 --simtype $3 -e 0.4 -z $4

echo "mean = 0.8, pi = 0.8, r=$1, nu=$2, $4 samples, 22 microbes per sample"
julia simulation/run_simulation.jl -n 2 -m 0.8 -p 0.8 -r $1 -k 22 -s 2358 -u $2 --simtype $3 -e 0.4 -z $4

sleep 4s

echo "mean = 1.6, pi = 0.3, r=$1, nu=$2, $4 samples, 8 microbes per sample"
julia simulation/run_simulation.jl -n 2 -m 1.6 -p 0.3 -r $1 -k 8 -s 2358 -u $2 --simtype $3 -e 0.4 -z $4

echo "mean = 1.6, pi = 0.3, r=$1, nu=$2, $4 samples, 22 microbes per sample"
julia simulation/run_simulation.jl -n 2 -m 1.6 -p 0.3 -r $1 -k 22 -s 2358 -u $2 --simtype $3 -e 0.4 -z $4

sleep 4s

echo "mean = 1.6, pi = 0.8, r=$1, nu=$2, $4 samples, 8 microbes per sample"
julia simulation/run_simulation.jl -n 2 -m 1.6 -p 0.8 -r $1 -k 8 -s 2358 -u $2 --simtype $3 -e 0.4 -z $4

echo "mean = 1.6, pi = 0.8, r=$1, nu=$2, $4 samples, 22 microbes per sample"
julia simulation/run_simulation.jl -n 2 -m 1.6 -p 0.8 -r $1 -k 22 -s 2358 -u $2 --simtype $3 -e 0.4 -z $4