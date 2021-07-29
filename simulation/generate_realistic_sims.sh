#!/bin/bash

echo "realistic simulations (sim 2)"

echo "Additive Phylogenetic"
bash simulation/generate_realistic_onetype.sh "additive_phylo"

echo "Additive Random"
bash simulation/generate_realistic_onetype.sh "additive_random"

echo "Interaction Phylogenetic"
bash simulation/generate_realistic_onetype.sh "interaction_phylo"

echo "Interaction Random"
bash simulation/generate_realistic_onetype.sh "interaction_random"

echo "Functional Redundancy Phylogenetic"
bash simulation/generate_realistic_onetype.sh "redundant_phylo"

echo "Functional Redundancy Random"
bash simulation/generate_realistic_onetype.sh "redundant_random"
