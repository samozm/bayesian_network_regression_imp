using ArgParse,Distributions,Random,CSV,RCall
include("../src/BayesNet.jl")

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
    "--nsamp", "-n"
        help = "number of samples to generate"
        arg_type = Int
        default = 70
    "--casenum", "-c"
        help = "case number to append to the name of the csv files create"
        arg_type = Int
        default = 1
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
    args_in = parse_CL_args()
    intxn = 0
    t = args_in["tottaxa"]
    k = args_in["samptaxa"]
    n = args_in["nsamp"]
    s = args_in["seed"]
    casen = args_in["casenum"]
    μₛ = args_in["mean"]
    πₛ = args_in["pi"]
    jcon = args_in["juliacon"]
    q = floor(Int,t*(t-1)/2)

    B,ξ = generate_Bs(t,μₛ=μₛ,πₛ=πₛ)

    y,A,m = generate_unrealistic_data(B,t,k,n,s)#,0,0.25)

    X = Matrix{Float64}(undef, size(A,1), q)
    for i in 1:size(A,1)
        X[i,:] = upper_triangle(A[i])
    end

    out_df = DataFrame(X,:auto)
    out_df[!,:y] = y

    if jcon
        CSV.write("juliacon/data/simulation1_case$(casen).csv",out_df)
        CSV.write("juliacon/data/simulation1_case$(casen)_bs.csv",DataFrame(B=upper_triangle(B)))
        CSV.write("juliacon/data/simulation1_case$(casen)_ms.csv",DataFrame(transpose(hcat(m...)),:auto))
        CSV.write("juliacon/data/simulation1_case$(casen)_xis.csv",DataFrame(TrueXi=ξ))
        return
    end
    CSV.write("data/simulation/simulation1_case$(casen).csv",out_df)
    CSV.write("data/simulation/simulation1_case$(casen)_bs.csv",DataFrame(B=upper_triangle(B)))
    CSV.write("data/simulation/simulation1_case$(casen)_ms.csv",DataFrame(transpose(hcat(m...)),:auto))
    CSV.write("data/simulation/simulation1_case$(casen)_xis.csv",DataFrame(TrueXi=ξ))

    #println("y")
    #show(stdout,"text/plain",y)
    #println("")
    #println("A[1]")
    #show(stdout,"text/plain",A[1])
    #println("")
    #println("m")
    #show(stdout,"text/plain",m)
    #println("")
end

function generate_Bs(t; μₛ=0.8,σₛ=1,πₛ=0.1)

    ξ = rand(Bernoulli(πₛ),t)#map(i -> rand(Normal(μₘ,σₘ)),1:t)
    B = zeros(t,t)

    for i in 1:t
        for j in (i+1):t
            if (ξ[i] == 1 && ξ[j] == 1)
                B[i,j] = B[j,i] = -abs(rand(Normal(μₛ,σₛ)))
            end
        end
    end

    return B,ξ
end



function generate_unrealistic_data(B,t,k,n,seed)#,μₐ,σₐ)

    Random.seed!(seed)
    y = zeros(n)
    m = [zeros(t)]
    A = [zeros(t,t)]
    A_base = zeros(t,t)
    @rput t
    R"source('src/sim_trees.R');inv_dist <- sim_tree_dists(t)"
    A_base = @rget inv_dist
    ϵ = rand(Normal(0,1),n)
    #for i in 1:t
        #for j in (i+1):t
        #    rnd = abs(rand(Normal(μₐ,σₐ)))
        #    while rnd == 0
        #        rnd = abs(rand(Normal(μₐ,σₐ)))
        #    end
        #    if rnd != 0
        #        rnd = 1/rnd #TODO: confirm we're gonna want to invert
        #    end
        #    A_base[i,j] = A_base[j,i] = rnd
        #end
    #end

    for i in 1:n
        chosen = sort(sample(1:t,k,replace=false))
        ind = zeros(n,n)
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
