#!/bin/bash

echo "mean = 0.8, pi = 0.3, samptaxa=8, type=$1, $2 samples"
julia simulation/realistic_sim.jl -s 734 -m 0.8 -p 0.3 -t 30 -k 8 --samplesize $2 --simtype $1

echo "mean = 0.8, pi = 0.3, samptaxa=22, type=$1, $2 samples"
julia simulation/realistic_sim.jl -s 734 -m 0.8 -p 0.3 -t 30 -k 22 --samplesize $2 --simtype $1 

echo "mean = 0.8, pi = 0.8, samptaxa=8, type=$1, $2 samples"
julia simulation/realistic_sim.jl -s 734 -m 0.8 -p 0.8 -t 30 -k 8 --samplesize $2 --simtype $1

echo "mean = 0.8, pi = 0.8, samptaxa=22, type=$1, $2 samples"
julia simulation/realistic_sim.jl -s 734 -m 0.8 -p 0.8 -t 30 -k 22 --samplesize $2 --simtype $1

echo "mean = 1.6, pi = 0.3, samptaxa=8, type=$1, $2 samples"
julia simulation/realistic_sim.jl -s 734 -m 1.6 -p 0.3 -t 30 -k 8 --samplesize $2 --simtype $1

echo "mean = 1.6, pi = 0.3, samptaxa=22, type=$1, $2 samples"
julia simulation/realistic_sim.jl -s 734 -m 1.6 -p 0.3 -t 30 -k 22 --samplesize $2 --simtype $1

echo "mean = 1.6, pi = 0.8, samptaxa=8, type=$1, $2 samples"
julia simulation/realistic_sim.jl -s 734 -m 1.6 -p 0.8 -t 30 -k 8 --samplesize $2 --simtype $1

echo "mean = 1.6, pi = 0.8, samptaxa=22, type=$1, $2 samples"
julia simulation/realistic_sim.jl -s 734 -m 1.6 -p 0.8 -t 30 -k 22 --samplesize $2 --simtype $1