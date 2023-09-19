include("realistic_sim.jl")
include("unrealistic_sim.jl")

g_rng = MersenneTwister(734)

@quickactivate

for πₑ=[0.3,0.8]
    for μₑ=[0.8,1.6]
        for n=[500,1000]
            for k=[8,22]
                for type=["additive_phylo","additive_random","interaction_phylo","interaction_random","redundant_phylo","redundant_random"]
                    rsd = sample(g_rng,100:9999)
                    generate_real(30,k,n,734,μₑ,πₑ,type,rsd)
                end
            end
        end
        for n=[100,500]
            for k=[8,15,22]
                rsd = sample(g_rng,100:9999)
                generate_unreal(30,k,n,734,μₑ,πₑ,false,rsd)
            end
        end
    end
end

for μₑ=[0.8,1.6]
    for πₑ=[0.0]
        for n=[100,500]
            for k=[8,15,22]
                rsd = sample(g_rng,100:9999)
                generate_unreal(30,k,n,734,μₑ,πₑ,false,rsd)
            end
        end
    end
end

for πₑ=[0.3,0.8]
    for μₑ=[0.8,1.6]
        for n=[50]
            for k=[8,22]
                for type=["additive_phylo","additive_random","interaction_phylo","interaction_random","redundant_phylo","redundant_random"]
                    rsd = sample(g_rng,100:9999)
                    generate_real(30,k,n,734,μₑ,πₑ,type,rsd)
                end
            end
        end
        for n=[50,500,1000]
            for k=[8,22,30]
                for type=["additive_phylo","additive_random","interaction_phylo","interaction_random","redundant_phylo","redundant_random"]
                    rsd = sample(g_rng,100:9999)
                    generate_real(60,k,n,734,μₑ,πₑ,type,rsd)
                end
            end
        end
end