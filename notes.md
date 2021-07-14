# Meeting with professor Rick Lankau

I met with Rick on April 6, 2021 to talk about potential ways to simulate data for microbiome.

He shared his review paper where he describes three 6 types of relationships between a "causal agent" x,z (microbial taxa) and a response y (phenotype). This is in Figure 1 of the BradfordHill paper.

The 6 types of relationships are:
- parallel (sub additive)
- additive
- synergistic (super additive)
- context-dependent
- indirect (mediation)
- spurious (confounding)

## Idea for simulation

1. Select t microbial taxa (say t=20)
2. Choose a phylogenetic tree that represents all n taxa
3. Simulate sequences on the phylogenetic tree
4. For each sample (i=1...n),
    - We will extract a random sample of k < t taxa
    - We will estimate the phylogenetic tree for the k taxa
    - We will calculate the distance matrix based on phylogenetic distances
    - We will create the microbiome network based on the distance matrix (if the adjacency matrix has to be binary, we can set a threshold to the distance so that if the distance < threshold, we set a 1). This is the network_i for the regression
    - We will set the response y as a function on the present taxa in the sample (see details below). This is the y_i for the regression
5. At the end, we have N_1,...N_n networks and y_1,...,y_n responses to run our Bayesian network regression model

### Response function
The response will be given by: y=b_1*m_1 + ... + b_k*m_k where m_1,...m_k are the microbial taxa present in the sample (0/1) and b_i is the effect of microbial taxon i.

We will do three causal relationship settings:
1. Additive model: y=b_1*m_1 + ... + b_k*m_k
2. Sub-additive model: y=b_1*m_1 + ... + b_k*m_k + Interaction where Interaction is negative
3. Super-additive model: y=b_1*m_1 + ... + b_k*m_k + Interaction where Interaction is positive

For each of these settings, we will choose the effects (b_i) and Interaction in two ways:
1. Randomly chosen numbers
2. Phylogenetically conserved: In this setup, we will simulate a trait value on the phylogenetic tree that will represent the effect of the microbes. The idea is that closely related microbes should have similar effect.


# Side interesting question
The effect of microbial taxa on phenotypes could be functionally redundant. That is, individually, each taxa can have an effect on the phenotype, but when put together, one of them stops having any effect on the phenotype.

One holy grail goal from machine-learning approaches could be to know which of the microbial taxa are functionally redundant from observational data (not experimental data).


# Follow-up meeting with Rick on simulations

Our process begins with a phylogenetic tree with 100 microbes, and then we take n samples from this tree with k microbes each. Each sample is one community.

The sampling process can be simple or complex depending on different factors:
- gradient (environmental filter): if all environments are expected to be similar, then all communities can be expected to be similar. If environments are expected to differ (in ph e.g.), then communities would be different and there would be biases on the microbes we observe on certain communities
- presence of other microbes: this affects mainly how the response is computed. For example, you have your pruned tree and you can check if it is overdispersed (species are less related than by chance: this is the case when you sample only one representative per clade) or underdispersed (species are more related than by chance: you have whole parts of the tree not included in the sample)
    - overdisperse case: species do not co-exist well, so you do not expect to see many interactions
    - underdisperse case: missing big chunks of the tree; trait filtering (you only observe microbes with certain trait) and interaction is very important

Sampling scheme:
1. Assume environments are similar, so we expect variation on communities only due to randon noise: it is ok to uniformly at random sample the microbes
2. Assume environments are different: in this case, some microbes will have a higher probability of appearing in certain environments

Computation of the response:
1. Phylogenetically informed microbe effects vs randomly chosen microbe effects
2. Additive model (ignore interactions) vs functional redundancy (there is a maximum value of the response) vs superadditive (some species together have a positive interaction)


## Things to keep in mind
- In real life, usually samples have a few common species very abundant, and many uncommon species (not as abundant)
- Note that we are not including any information on abundances on samples (because we do not include this in the model). This could be an extension of the model