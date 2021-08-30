using ArgParse,Distributions,Random,CSV,RCall,DrWatson
using LinearAlgebra,BayesianNetworkRegression,DataFrames
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
        default = 100
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
    μₛ = args_in["mean"]
    πₛ = args_in["pi"]
    jcon = args_in["juliacon"]
    q = floor(Int,t*(t-1)/2)

    Random.seed!(seed)
    B,ξ = generate_Bs(t,μₛ=μₛ,πₛ=πₛ)

    y,A,m = generate_unrealistic_data(B,t,k,n,seed)#,0,0.25)

    X = Matrix{Float64}(undef, size(A,1), q)
    for i in 1:size(A,1)
        X[i,:] = BayesianNetworkRegression.lower_triangle(A[i])
    end

    out_df = DataFrame(X,:auto)
    out_df[!,:y] = y

    saveinfo = Dict("simnum"=>"1","pi"=>πₛ,"mu"=>μₛ,"n_microbes"=>k,"samplesize"=>n)
    output_data(saveinfo,out_df,B,m,ξ,jcon,"unrealistic")

end

function generate_Bs(t; μₛ=0.8,σₛ=1,πₛ=0.1)

    ξ = rand(Bernoulli(πₛ),t)
    B = zeros(t,t)

    fill_B(B,ξ,t,μₛ,σₛ)

    return B,ξ
end

function generate_unrealistic_data(B,t,k,n,seed)
    y = zeros(n)
    m = [zeros(t)]
    A = [zeros(t,t)]
    A_base = zeros(t,t)
    @rput t
    @rput seed
    R"set.seed(seed);source('src/sim_trees.R');inv_dist <- sim_tree_dists(t)"
    A_base = @rget inv_dist
    ϵ = rand(Normal(0,1),n)

    for i in 1:n
        chosen = sort(sample(1:t,k,replace=false))
        for j in chosen
            for l in chosen
                A[i][j,l] = A[i][l,j] = A_base[j,l]
            end
        end
        m[i][chosen] .= 1
        y[i] = tr(transpose(B) * A[i]) + ϵ[i]
        if i != n
            append!(A,[zeros(t,t)])
            append!(m,[zeros(t)])
        end
    end


    return y,A,m,ϵ
end

main()
