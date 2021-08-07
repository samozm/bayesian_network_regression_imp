using Base: String, Int64, simd_outer_range
using ArgParse,Distributions,Random,CSV,RCall,DrWatson
using LinearAlgebra,BayesianNetworkRegression,DataFrames,PhyloNetworks
include("sim_utils.jl")

function parse_CL_args()
    args = ArgParseSettings()
    @add_arg_table! args begin
    "--seed", "-s"
        help = "random seed to use for simulations"
        arg_type = Int
        default = 123
    "--tottaxa", "-t"
        help = "number of microbial taxa to generate"
        arg_type = Int
        default = 100
    "--samptaxa","-k"
        help = "number of taxa to actually use in each sample"
        arg_type = Int
        default = 20
    "--samplesize", "-z"
        help = "number of samples to generate"
        arg_type = Int
        default = 70
    "--mean", "-m"
        help = "mean to use to draw edge coefficients"
        arg_type = Float64
        default = 0.8
    "--pi", "-p"
        help = "probability of any single node being influential on the response"
        arg_type = Float64
        default = 0.1
    "--juliacon", "-j"
        help = "flag indicating output should go to juliacon folder"
        action = :store_true
    "--simtype", "-y"
        help = "type of simulation to run: additive_phylo, additive_random, interaction_phylo, interaction_random, redundant_phylo, or redundant_random"
        arg_type = String
        default = "additive_phylo"
    end
    return parse_args(args)
end


function main()
    @quickactivate
    args_in = parse_CL_args()
    t = args_in["tottaxa"]
    k = args_in["samptaxa"]
    n = args_in["samplesize"]
    seed = args_in["seed"]
    μₑ = args_in["mean"]
    πₑ = args_in["pi"]
    jcon = args_in["juliacon"]
    type = args_in["simtype"]
    q = floor(Int,t*(t-1)/2)

    Random.seed!(seed)

    μₛ = 0.4
    σₛ = 1.0

    ξ,B,y,A,m = generate_realistic_data(t,k,n,μₑ,πₑ,μₛ,σₛ,type)#,0,0.25)

    X = Matrix{Float64}(undef, size(A,1), q)
    for i in 1:size(A,1)
        X[i,:] = BayesianNetworkRegression.lower_triangle(A[i])
    end

    out_df = DataFrame(X,:auto)
    out_df[!,:y] = y

    saveinfo = Dict("simnum"=>"2","pi"=>πₑ,"mu"=>μₑ,"n_microbes"=>k,"type"=>type,"edge_mu"=>μₛ,"samplesize"=>n)
    output_data(saveinfo,out_df,B,m,ξ,jcon,"realistic")

    saveinfo["out"] = "main-effects"
    CSV.write(datadir(joinpath("simulation","realistic"),savename(saveinfo,"csv",digits=1)),DataFrame(me=diag(B)))

end

function generate_realistic_data(t,k,n,μₑ,πₑ,μₛ,σₛ,type)

    @rput t
    R"source('src/sim_trees.R');tree <- sim_tree_string(t); A <- tree_dist(tree$tree); tree_str <- tree$tree_str"
    tree = @rget tree_str
    A_base = @rget A

    ξ = rand(Bernoulli(πₑ),t)

    type_arr = split(type,"_")

    L_dict = Dict(0.8 => 1, 1.6 => 2)

    if type_arr[2] == "phylo"
        # generate C (main effects)
        tree_obj = readTopology(tree)
        trait_params = ParamsBM(μₛ,σₛ)
        tree_sim = simulate(tree_obj, trait_params)
        C = tree_sim[:Tips]
        B = diagm(C)

        if type_arr[1] == "additive"
            return generate_yAm(ξ,B,t,k,n,A_base)
        end

        fill_B(B,ξ,t,μₛ,σₛ)

        if type_arr[1] == "interaction"
            return generate_yAm(ξ,B,t,k,n,A_base)
        elseif type_arr[1] == "redundant" 
            return generate_yAm(ξ,B,t,k,n,A_base,true,L=L_dict[μₑ])
        end

    elseif type_arr[2] == "random"
        # generate C (main effects)
        C = rand(Normal(μₑ,1),t)
        B = diagm(C)

        if type_arr[1] == "additive"
            return generate_yAm(ξ,B,t,k,n,A_base)
        end

        fill_B(B,ξ,t,μₛ,σₛ)

        if type_arr[1] == "interaction" 
            return generate_yAm(ξ,B,t,k,n,A_base)
        elseif type_arr[1] == "redundant" 
            return generate_yAm(ξ,B,t,k,n,A_base,true,L=L_dict[μₑ])
        end
    end
end


function generate_yAm(ξ,B,t,k,n,A_base,redundant=false;L=1)
    y = zeros(n)
    m = [zeros(t)]
    A = [zeros(t,t)]
    ϵ = rand(Normal(0,1),n)


    for i in 1:n
        chosen = sort(sample(1:t,k,replace=false))
        for j in chosen
            for l in chosen
                A[i][j,l] = A[i][l,j] = A_base[j,l]
            end
        end
        m[i][chosen] .= 1
        #TODO - make sure this is right (IND[i,j] should be 1 if m[i] and ξ[j] are 1, otherwise 0)
        IND = m[i] * transpose(ξ)
        if redundant
            y[i] = max(tr(transpose(B) * IND) + ϵ[i],L)
        else
            y[i] = tr(transpose(B) * IND) + ϵ[i]
        end
        if i != n
            append!(A,[zeros(t,t)])
            append!(m,[zeros(t)])
        end
    end

    return ξ,B,y,A,m
end

main()