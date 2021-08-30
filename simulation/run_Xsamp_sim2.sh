#!/bin/bash

echo "realistic simulations (sim 2) - $1 samples"

#additive_phylo, additive_random, interaction_phylo, interaction_random, redundant_phylo, or redundant_random
echo "############# additive_phylo #############"
bash simulation/run_realistic.sh 9 10 "additive_phylo" $1

sleep 4s

echo "############# additive_random #############"
bash simulation/run_realistic.sh 9 10 "additive_random" $1

sleep 4s

echo "############# interaction_phylo #############"
bash simulation/run_realistic.sh 9 10 "interaction_phylo" $1

sleep 4s

echo "############# interaction_random #############"
bash simulation/run_realistic.sh 9 10 "interaction_random" $1

sleep 4s

echo "############# redundant_phylo #############"
bash simulation/run_realistic.sh 9 10 "redundant_phylo" $1

sleep 4s

echo "############# redundant_random #############"
bash simulation/run_realistic.sh 9 10 "redundant_random" $1