#!/bin/bash

echo "realistic simulations (sim 2)"

echo "Additive Phylogenetic"
bash simulation/generate_realistic_onetype.sh "additive_phylo" $1

echo "Additive Random"
bash simulation/generate_realistic_onetype.sh "additive_random" $1

echo "Interaction Phylogenetic"
bash simulation/generate_realistic_onetype.sh "interaction_phylo" $1

echo "Interaction Random"
bash simulation/generate_realistic_onetype.sh "interaction_random" $1

echo "Functional Redundancy Phylogenetic"
bash simulation/generate_realistic_onetype.sh "redundant_phylo" $1

echo "Functional Redundancy Random"
bash simulation/generate_realistic_onetype.sh "redundant_random" $1
