include ArgParse,Distributions,Random

function parse_CL_args()
    args = ArgParseSettings()
    @add_arg_table! args begin
    "--seed", "-s"
        help = "random seed to use for simulations"
        argtype = Int
    "-t"
        help = "number of microbial taxa to use"
        argtype = Int
        default = 20
    "--nsamp", "-n"
        help = "number of samples to generate"
        argtype = Int
        default = 70
    end
    return parse_args(args)
end


function main()

end


function generate_unrealistic_data(B,t,k,n,seed,interaction)

    chosen_m = zeros(n)
    y = zeros(n)
    m = map(i -> sample([0,1]),1:t)

    if(size(B)!= t && interaction==0)
        println("Size of B does not match number of taxa")
        return
    end

    if interaction == -1
        for i in 1:n
            chosen_m[i] = sample(m,k,replace=false)
            y[i] = chosen_m * transpose(B)
        end
    elseif interaction == 1
        for i in 1:n
            chosen_m[i] = sample(m,k,replace=false)
            y[i] = chosen_m * transpose(B)
        end

    elseif interaction == 0
        for i in 1:n
            chosen_m[i] = sample(m,k,replace=false)
            y[i] = chosen_m * transpose(B)
        end
    end
    return y,chosen_m
end
