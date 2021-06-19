using CSV,ArgParse,TickTock
include("../src/BayesNet.jl")
include("../src/plot_output.jl")


function parse_CL_args()
    args = ArgParseSettings()
    @add_arg_table! args begin
    "--nburn", "-b"
        help="Number of burn-in Gibbs samples to take"
        arg_type = Int
        default = 30000
    "--nsamp", "-s"
        help="Number of Gibbs samples to keep (after burn-in)"
        arg_type = Int
        default = 20000
    "--simnum", "-n"
        help="which simulation to run"
        arg_type = Int
        default = 1
    "--casenum", "-c"
        help="which case of the simulation to run"
        arg_type = Int
        default = 1
    end
    return parse_args(args)
end

function main()
    parsed_CL_args = parse_CL_args()
    nburn = parsed_CL_args["nburn"]
    nsamp = parsed_CL_args["nsamp"]
    simnum = parsed_CL_args["simnum"]
    casenum = parsed_CL_args["casenum"]
    γ,MSE,ξ = sim_one_case(simnum,casenum,nburn,nsamp)

    output_results(γ[:,nburn+1:nburn+nsamp],MSE,mean(ξ[nburn+1:nburn+nsamp]),simnum,casenum)
end

function sim_one_case(simnum,casenum,nburn,nsamp,η=1.01,ζ=1,ι=1,R=5,aΔ=1,bΔ=1,ν=10)
    data_in = DataFrame(CSV.File("data/simulation/simulation$(simnum)_case$(casenum).csv"))

    X = convert(Matrix,data_in[!,names(data_in,Not("y"))])
    y = data_in[:,:y]

    b_in = DataFrame(CSV.File("data/simulation/simulation$(simnum)_case$(casenum)_bs.csv"))
    B₀ = convert(Array{Float64,1},b_in[!,:B])

    q = size(X,2)
    V = convert(Int,(1 + sqrt(1 + 8*q))/2)

    tick()
    τ², u, ξ, γ, D, θ, Δ, M, μ, Λ, πᵥ = BayesNet(X, y, R, nburn=nburn,nsamples=nsamp, V_in = V, aΔ=aΔ, bΔ=bΔ,ν=ν,ι=ι,ζ=ζ,x_transform=false)
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

function output_results(γ,MSE,ξ,simnum,casenum)
    q = size(γ,1)
    V = convert(Int,(1 + sqrt(1 + 8*q))/2)

    plot_γ_sim(γ, "Gamma", "simulation$(simnum)_case$(casenum)")
    show(stdout,"text/plain",DataFrame(Xi=ξ))
    println("")
    show(stdout,"text/plain",DataFrame(MSE=MSE))
    println("")
    println("Gamma")
    gam = DataFrame(create_upper_tri(vec(median.(γ)),V),:auto)
    show(stdout,"text/plain",gam[!,names(gam,Not("x1"))])
    println("")
end

main()
