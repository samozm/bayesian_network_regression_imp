
echo "realistic simulations (sim 2)"

echo "Additive Phylogenetic"
bash simulation/generate_realistic_onetype.sh "additive_phylo"

echo "Additive Random"
bash simulation/generate_realistic_onetype.sh "additive_random"

echo "Interaction Phylogenetic"
bash simulation/generate_realistic_onetype.sh "interaction_phylo"

echo "Interaction Random"
bash simulation/generate_realistic_onetype.sh "interaction_random"

echo "Functional Redundancy Phylogenetic L=1"
bash simulation/generate_realistic_onetype.sh "redundant_phylo" 1

echo "Functional Redundancy Random L=1"
bash simulation/generate_realistic_onetype.sh "redundant_random" 1

echo "Functional Redundancy Phylogenetic L=2"
bash simulation/generate_realistic_onetype.sh "redundant_phylo" 2

echo "Functional Redundancy Random L=2"
bash simulation/generate_realistic_onetype.sh "redundant_random" 2
