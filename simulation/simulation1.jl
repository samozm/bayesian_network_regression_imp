using ArgParse,Distributions,Random,CSV
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
        default = 20
    "--samptaxa","-k"
        help = "number of taxa to actually use in each sample"
        arg_type = Int
        default = 10
    "--nsamp", "-n"
        help = "number of samples to generate"
        arg_type = Int
        default = 70
    "--interaction", "-i"
        help = "interaction effect direction between microbes: 0,1,-1"
        arg_type = Int
        default = 0
    "--simnum", "-m"
        help = "number to append to the name of the csv files create"
        arg_type = Int
        default = 1
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
    simn = args_in["simnum"]
    q = floor(Int,t*(t-1)/2)
    if !(args_in["interaction"] in [0,1,-1])
        println("invalid interaction value given. using no interaction (0)")
    else
        intxn = args_in["interaction"]
    end

    B,M_eff = generate_Bs(intxn,t)

    y,A,m = generate_unrealistic_data(B,M_eff,t,k,n,s,0,0.25)

    X = Matrix{Float64}(undef, size(A,1), q)
    for i in 1:size(A,1)
        X[i,:] = upper_triangle(A[i])
    end

    out_df = DataFrame(X,:auto)
    out_df[!,:y] = y

    CSV.write("data/simulation/simulation$(simn).csv",out_df)
    CSV.write("data/simulation/simulation$(simn)_bs.csv",DataFrame(B=upper_triangle(B)))
    CSV.write("data/simulation/simulation$(simn)_ms.csv",DataFrame(transpose(hcat(m...)),:auto))

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

function generate_Bs(intxn,t,μₘ=0,σₘ=1,μ_b=0.5,σ_b=1,p=0.8)

    M_eff = zeros(t)#map(i -> rand(Normal(μₘ,σₘ)),1:t)
    B = zeros(t,t)

    if intxn == -1
        for i in 1:t
            if rand(Bernoulli(p))
                M_eff[i] = rand(Normal(μₘ,σₘ))
            end
            for j in (i+1):t
                B[i,j] = B[j,i] = -abs(rand(Normal(μ_b,σ_b)))
            end
        end
    elseif intxn == 1
        for i in 1:t
            for j in (i+1):t
                B[i,j] = B[j,i] = abs(rand(Normal(μ_b,σ_b)))
            end
        end
    end

    return B,M_eff
end



function generate_unrealistic_data(B,M_eff,t,k,n,seed,μₐ,σₐ)

    Random.seed!(seed)
    y = zeros(n)
    m = [zeros(t)]
    A = [zeros(t,t)]
    A_base = zeros(t,t)
    for i in 1:t
        for j in (i+1):t
            rnd = abs(rand(Normal(μₐ,σₐ)))
            while rnd == 0
                rnd = abs(rand(Normal(μₐ,σₐ)))
            end
            if rnd != 0
                rnd = 1/rnd #TODO: confirm we're gonna want to invert
            end
            A_base[i,j] = A_base[j,i] = rnd
        end
    end

    for i in 1:n
        chosen = sort(sample(1:t,k,replace=false))
        ind = zeros(n,n)
        for j in 1:t
            for l in (j+1):t
                if j in chosen && l in chosen
                    A[i][j,l] = A[i][l,j] = A_base[j,l]
                end
            end
        end
        m[i][chosen] .= 1
        y[i] = transpose(m[i]) * M_eff + tr(transpose(B) * A[i])
        if i != n
            append!(A,[zeros(t,t)])
            append!(m,[zeros(t)])
        end
    end


    return y,A,m
end

main()
