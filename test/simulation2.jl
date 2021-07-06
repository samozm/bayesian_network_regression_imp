using Distributions,Random,TickTock,CSV,ArgParse,Printf,RCall
include("../src/BayesNet.jl")
include("../src/plot_output.jl")


function parse_CL_args()
    args = ArgParseSettings()
    @add_arg_table! args begin
        "--case", "-c"
            help="Case number, 1 or 2 controls sparsity parameter"
            arg_type = Int
            default = 1
        "--nburn", "-b"
            help="Number of burn-in Gibbs samples to take"
            arg_type = Int
            default = 30000
        "--nsamp", "-s"
            help="Number of Gibbs samples to keep (after burn-in)"
            arg_type = Int
            default = 20000
    end

    return parse_args(args)
end

function main()
    parsed_CL_args = parse_CL_args()
    case = parsed_CL_args["case"]
    nburn = parsed_CL_args["nburn"]
    nsamp = parsed_CL_args["nsamp"]
    γ₁,MSE₁,ξ₁ = sim_one_case(case,nburn,nsamp)
    @rput nburn
    @rput nsamp
    if (case == 1)
        R"source('test/run_guha_sim2_case1.R')"
        println("Guha:")
        tick()
        R"retlist <- guha_sim2_case1(nburn,nsamp)"
        tock()
    else
        R"source('test/run_guha_sim2_case2.R')"
        println("Guha:")
        tick()
        R"retlist <- guha_sim2_case2(nburn,nsamp)"
        tock()
    end
    R"gamma <- retlist$gamma; MSE <- retlist$MSE; xis <- retlist$xis"
    γ₂ = @rget gamma
    MSE₂ = @rget MSE
    ξ₂ = @rget xis
    output_results(γ₁[:,nburn+1:nburn+nsamp],MSE₁,mean(ξ₁[nburn+1:nburn+nsamp]),γ₂,MSE₂,ξ₂,case)
end

function output_results(γ₁,MSE₁,ξ₁,γ₂,MSE₂,ξ₂,case)
    plot_γs_test(γ₁, γ₂, "Mine", "Guha","sim2_case$case")
    show(stdout,"text/plain",DataFrame(My_Xi=ξ₁,Guha_Xi=ξ₂))
    println("")
    show(stdout,"text/plain",DataFrame(My_MSE=MSE₁,Guha_MSE=MSE₂))
    println("")
end

function sim_one_case(case,nburn,nsamp)
    η  = 1.01
    ζ  = 1.0
    ι  = 1.0
    R  = 5
    aΔ = 1.0
    bΔ = 1.0
    ν = 10
    V = 20
    q = floor(Int,V*(V-1)/2)

    data_in = DataFrame(CSV.File("data/test/simulation2_case$(case).csv"))

    X = Matrix(data_in[:,1:190])
    y = data_in[:,191]

    b_in = DataFrame(CSV.File("data/test/simulation2_case$(case)_bs.csv"))
    B₀ = convert(Array{Float64,1},b_in[!,:B])

    @printf("Sim 2 Case %d",case)
    println("")

    tick()
    τ², u, ξ, γ, D, θ, Δ, M, μ, Λ, πᵥ = BayesNet(X, y, R, nburn=nburn,nsamples=nsamp, V_in = 20, aΔ=1.0, bΔ=1.0,ν=10,ι=1.0,ζ=1.0,x_transform=false)
    tock()

    low = zeros(190)
    high = zeros(190)
    lw = convert(Int64, round(nsamp * 0.1))
    hi = convert(Int64, round(nsamp * 0.9))

    γ_n = hcat(γ...)

    γ_n2 = mean(γ[nburn+1:nburn+nsamp])
    γ₀ = B₀
    MSE = 0
    for i in 1:190
        MSE = MSE + (γ_n2[i] - γ₀[i])^2
    end
    MSE = MSE * (2/(V*(V-1)))

    return γ_n,MSE,ξ
end


main()
