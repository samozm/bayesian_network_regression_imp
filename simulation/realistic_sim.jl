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
