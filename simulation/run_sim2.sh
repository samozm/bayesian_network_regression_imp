#!/bin/bash

echo "realistic simulations (sim 2)"

sleep 3600s

#additive_phylo, additive_random, interaction_phylo, interaction_random, redundant_phylo, or redundant_random
echo "############# additive_phylo #############"
bash simulation/run_unrealistic.sh 9 10 "additive_phylo"

sleep 4s

echo "############# additive_random #############"
bash simulation/run_unrealistic.sh 9 10 "additive_random"

sleep 4s

echo "############# interaction_phylo #############"
bash simulation/run_unrealistic.sh 9 10 "interaction_phylo"

sleep 4s

echo "############# interaction_random #############"
bash simulation/run_unrealistic.sh 9 10 "interaction_random"

sleep 4s

echo "############# redundant_phylo #############"
bash simulation/run_unrealistic.sh 9 10 "redundant_phylo"

sleep 4s

echo "############# redundant_random #############"
bash simulation/run_unrealistic.sh 9 10 "redundant_random"