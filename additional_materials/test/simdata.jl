using CSV,Random,Distributions,LinearAlgebra,Printf
include("../src/BayesNet.jl")


function parse_CL_args()
    args = ArgParseSettings()
    @add_arg_table! args begin
        "--case", "-c"
            help="Case number, 1 or 2 controls sparsity parameter"
            arg_type = Int
            default = 1
    end

    return parse_args(args)
end

function main()
    parsed_CL_args = parse_CL_args()
    generate_data(parsed_CL_args["case"])
end

function generate_data(case)

    η  = 1.01
    ζ  = 1
    ι  = 1
    R  = 5
    aΔ = 1
    bΔ = 1
    ν = 10
    V = 20
    q = floor(Int,V*(V-1)/2)

    n=70
    π₂=0.3

    if π₂==0.8
        println("Case 2")
    elseif π₂==0.3
        println("Case 1")
    else
        println("Unknown pi value")
    end
    A  = [zeros(V,V)]
    y = zeros(n)

    ξ⁰ = map(l -> rand(Bernoulli(π₂)),1:V)
    B₀ = zeros(V,V)

    for k = 1:V
        for l = (k+1):V
            if ξ⁰[k] == 1 && ξ⁰[l] == 1
                B₀[k,l] = rand(Normal(0.8,1))
                B₀[l,k] = B₀[k,l]
            else
                B₀[k,l] = 0
                B₀[l,k] = 0
            end
        end
    end

    for i = 1:n
        for k=1:V
            for l = (k+1):V
                A[i][k,l] = rand(Normal(0,1))
                A[i][l,k] = A[i][k,l]
            end
        end
        τ₀² = 1
        ϵᵢ = rand(Normal(0,τ₀²))
        y[i] = tr(transpose(B₀) * A[i]) + ϵᵢ
        #y[i] = sum(upper_triangle(B₀) .* upper_triangle(A[i])) + ϵᵢ
        if i != n
            append!(A,[zeros(V,V)])
        end
    end

    X_new = Matrix{Float64}(undef, size(A,1), q)
    for i in 1:size(A,1)
        X_new[i,:] = upper_triangle(A[i])
    end


    out_df = DataFrame(X_new,:auto)
    out_df[!,:y] = y

    CSV.write("data/test/simulation2_case$(case).csv",out_df)
    CSV.write("data/test/simulation2_case$(case)_bs.csv",DataFrame(B=upper_triangle(B₀)))
end
