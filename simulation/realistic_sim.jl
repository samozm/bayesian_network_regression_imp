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
    "--gseed", "-g"
        help = "random seed used to generate gaussian noise. different from seed so that the same seed can be used between simulation types but get different gaussian noise"
        arg_type = Int
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
    type = args_in["simtype"]
    gseed = args_in["gseed"]

    generate_real(t,k,n,seed,μₑ,πₑ,type,gseed)
end

function generate_real(t,k,n,seed,μₑ,πₑ,type,gseed,normdiag=false)
    q = floor(Int,t*(t+1)/2)

    Random.seed!(seed)
    rng = MersenneTwister(seed)
    gauss_rng = MersenneTwister(gseed)

    μₛ = 0.4
    σₛ = 1.0

    aug = false
    num_aug=0

<<<<<<< HEAD
    if n == 100
=======
    if n == 100 
>>>>>>> d5f23a3283dc0459a2e4ead6e8e867926119458d
        old_n = n
        n = 50
        aug = true
        num_aug=1
        loadinfo = Dict("simnum"=>"2","pi"=>πₑ,"mu"=>μₑ,"n_microbes"=>k,"type"=>type,"edge_mu"=>μₛ,"samplesize"=>n)
        loadinfo["out"] = "XYs"
        data_in = DataFrame(CSV.File(string("data/simulation/realistic/",savename(loadinfo,"csv",digits=1))))
        X = Matrix(data_in[:,names(data_in,Not("y"))])
        y = Vector(data_in[:,:y])

        loadinfo["out"] = "xis"
        ξ_in = DataFrame(CSV.File(string("data/simulation/realistic/",savename(loadinfo,"csv",digits=1))))[!,"TrueXi"]

        loadinfo["out"] = "A"
        A_base = DataFrame(CSV.File(string("data/simulation/realistic/",savename(loadinfo,"csv",digits=1))))
        A_base = Matrix(A_base)

        X,y = augment(X,y,ξ_in,A_base,rng,num_aug)
        n = old_n
        ξ = ξ_in

        loadinfo["out"] = "bs"
        B = DataFrame(CSV.File(string("data/simulation/realistic/",savename(loadinfo,"csv",digits=1))))

    elseif n == 200
        old_n = n
        n = 50
        aug = true
        num_aug=3
        loadinfo = Dict("simnum"=>"2","pi"=>πₑ,"mu"=>μₑ,"n_microbes"=>k,"type"=>type,"edge_mu"=>μₛ,"samplesize"=>n)
        loadinfo["out"] = "XYs"
        data_in = DataFrame(CSV.File(string("data/simulation/realistic/",savename(loadinfo,"csv",digits=1))))
        X = Matrix(data_in[:,names(data_in,Not("y"))])
        y = Vector(data_in[:,:y])

        loadinfo["out"] = "xis"
        ξ_in = DataFrame(CSV.File(string("data/simulation/realistic/",savename(loadinfo,"csv",digits=1))))[!,"TrueXi"]

        loadinfo["out"] = "A"
        A_base = DataFrame(CSV.File(string("data/simulation/realistic/",savename(loadinfo,"csv",digits=1))))
        A_base = Matrix(A_base)

        X,y = augment(X,y,ξ_in,A_base,rng,num_aug)
        n = old_n
        ξ = ξ_in

        loadinfo["out"] = "bs"
        B = DataFrame(CSV.File(string("data/simulation/realistic/",savename(loadinfo,"csv",digits=1))))
    else
        
        ξ,B,y,A,m,A_base = generate_realistic_data(t,k,n,μₑ,πₑ,μₛ,σₛ,type,seed,rng,gauss_rng,normdiag)#,0,0.25)

        X = Matrix{Float64}(undef, size(A,1), q)
        for i in 1:size(A,1)
            X[i,:] = lower_triangle(A[i])
        end
        B = DataFrame(B=lower_triangle(B)[:,1])
    end

    saveinfo = Dict("simnum"=>"2","pi"=>πₑ,"mu"=>μₑ,"n_microbes"=>k,"type"=>type,"edge_mu"=>μₛ,"samplesize"=>n)
    out_df = DataFrame(X,:auto)
    out_df[!,:y] = y

    output_data(saveinfo,out_df,ξ,A_base,B,false,"realistic")

end

function generate_realistic_data(t,k,n,μₑ,πₑ,μₛ,σₛ,type,seed,rng,gauss_rng,normdiag=false)

    @rput t
    @rput seed
    R"set.seed(seed);source('simulation/sim_trees.R');tree <- sim_tree_string(t); A <- tree_dist(tree$tree); tree_str <- tree$tree_str"
    tree = @rget tree_str
    A_base = @rget A
    mean_A_base = mean(A_base)

    if normdiag
        for i in 1:size(A_base,1)
            A_base[i,i] = mean_A_base
        end
    else    
        for i in 1:size(A_base,1)
            A_base[i,i] = 1
        end
    end


    ξ = rand(rng,Bernoulli(πₑ),t)

    type_arr = split(type,"_")

    # μₑ + πₑ
    L_dict = Dict(1.1 => 3,  1.6 => 22, 1.9 => 7, 2.4 => 30)

    if type_arr[2] == "phylo"
        # generate C (main effects)
        tree_obj = readTopology(tree)
        trait_params = ParamsBM(μₑ,σₛ)
        tree_sim = simulate(tree_obj, trait_params)
        C = tree_sim[:Tips]
        B = diagm(C) .* ξ

        if type_arr[1] == "additive"
            return generate_yAm(ξ,B,t,k,n,A_base,rng,gauss_rng)
        end

        fill_B(B,ξ,t,μₛ,σₛ,rng)

        if type_arr[1] == "interaction"
            return generate_yAm(ξ,B,t,k,n,A_base,rng,gauss_rng)
        elseif type_arr[1] == "redundant" 
            return generate_yAm(ξ,B,t,k,n,A_base,rng,gauss_rng,true,L=L_dict[round(μₑ+πₑ,digits=1)])
        end

    elseif type_arr[2] == "random"
        # generate C (main effects)
        C = rand(rng,Normal(μₑ,1),t)
        B = diagm(C) .* ξ

        if type_arr[1] == "additive"
            return generate_yAm(ξ,B,t,k,n,A_base,rng,gauss_rng)
        end

        fill_B(B,ξ,t,μₛ,σₛ,rng)

        if type_arr[1] == "interaction" 
            return generate_yAm(ξ,B,t,k,n,A_base,rng,gauss_rng)
        elseif type_arr[1] == "redundant" 
            return generate_yAm(ξ,B,t,k,n,A_base,rng,gauss_rng,true,L=L_dict[round(μₑ+πₑ,digits=1)])
        end
    end
end


function generate_yAm(ξ,B,t,k,n,A_base,rng,gauss_rng,redundant=false;L=1)
    y = zeros(n)
    m = [zeros(t)]
    A = [zeros(t,t)]
    #ϵ = rand(rng,Normal(0,1),n)


    for i in 1:n
        chosen = sort(sample(rng,1:t,k,replace=false))
        for j in chosen
            for l in chosen
                A[i][j,l] = A[i][l,j] = A_base[j,l]
            end
        end
        m[i][chosen] .= 1
       
        if redundant
            y[i] = min(tr(transpose(B) * A[i]) + rand(gauss_rng,Normal(0,1)),L)
        else
            y[i] = tr(transpose(B) * A[i]) + rand(gauss_rng,Normal(0,1))
        end
        if i != n
            append!(A,[zeros(t,t)])
            append!(m,[zeros(t)])
        end
    end

    return ξ,B,y,A,m,A_base
end


#main()