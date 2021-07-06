using LinearAlgebra: Matrix
using Base: Bool, Float16, Int16
using CSV,ArgParse,TickTock
using ProfileView, Traceur
using BayesianNetworkRegression
using Random, DataFrames, LinearAlgebra, StatsBase
using Distributions
using BenchmarkTools,ProfileView
#include("../src/BayesNet.jl")
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
    "--samptaxa","-k"
        help = "number of taxa actually used in each sample (relates to filename)"
        arg_type = Int
        default = 20
    "--mean", "-m"
        help = "mean to use to draw edge coefficients (relates to filename)"
        arg_type = Float64
        default = 0.8
    "--pi", "-p"
        help = "probability of any single node being influential on the response (relates to filename)"
        arg_type = Float64
        default = 0.1
    "-r"
        help = "r value, the dimension of the latent space (relates to filename)"
        arg_type = Int
        default = 5
    "--juliacon", "-j"
        help = "flag indicating output should go to juliacon folder"
        action = :store_true
    end
    return parse_args(args)
end

function main()
    parsed_CL_args = parse_CL_args()
    nburn = parsed_CL_args["nburn"]
    nsamp = parsed_CL_args["nsamp"]
    simnum = parsed_CL_args["simnum"]
    jcon = parsed_CL_args["juliacon"]
    μₛ = parsed_CL_args["mean"]
    πₛ = parsed_CL_args["pi"]
    R = parsed_CL_args["r"]
    k = parsed_CL_args["samptaxa"]

    run_case_and_output(nburn,nsamp,simnum,μₛ,πₛ,R,k,jcon)
end

function run_case_and_output(nburn,nsamp,simnum,μₛ,πₛ,R,k,jcon)
    γ,MSE,ξ = sim_one_case(simnum,nburn,nsamp,μₛ,πₛ,k,jcon,R=R)

    if jcon
        ξ_in = DataFrame(CSV.File("juliacon/data/sim1_pi$(πₛ)-mu$(μₛ)-$(k)microbes_xis.csv"))
    else
        ξ_in = DataFrame(CSV.File("data/simulation/sim1_pi$(πₛ)-mu$(μₛ)-$(k)microbes_xis.csv"))

        step = Int(floor(nsamp/100))
        if step==0
            step=1
        end
    end

    output_results(γ[1:nsamp],MSE,mean(ξ[1:nsamp]),ξ_in,simnum,πₛ,μₛ,R,k,jcon)
end

function sim_one_case(simnum,nburn,nsamp,μₛ,πₛ,k,jcon::Bool;η=1.01,ζ=1.0,ι=1.0,R=5,aΔ=1.0,bΔ=1.0,ν=10)
    if jcon
        data_in = DataFrame(CSV.File("juliacon/data/sim$(simnum)_pi$(πₛ)-mu$(μₛ)-$(k)microbes_XYs.csv"))

        b_in = DataFrame(CSV.File("juliacon/data/sim$(simnum)_pi$(πₛ)-mu$(μₛ)-$(k)microbes_bs.csv"))
        B₀ = convert(Array{Float64,1},b_in[!,:B])
    else
        data_in = DataFrame(CSV.File("data/simulation/sim$(simnum)_pi$(πₛ)-mu$(μₛ)-$(k)microbes_XYs.csv"))

        b_in = DataFrame(CSV.File("data/simulation/sim$(simnum)_pi$(πₛ)-mu$(μₛ)-$(k)microbes_bs.csv"))
        B₀ = convert(Array{Float64,1},b_in[!,:B])
    end
    #X = convert(Matrix,data_in[:,names(data_in,Not("y"))])
    X = Matrix(data_in[:,names(data_in,Not("y"))])
    y = data_in[:,:y]

    q = size(X,2)
    V = convert(Int,(1 + sqrt(1 + 8*q))/2)

    tick()
    #τ², u, ξ, γ, D, θ, Δ, M, μ, Λ, πᵥ = BayesNet(X, y, R, η=η, nburn=nburn,nsamples=nsamp, V_in = V, aΔ=aΔ, bΔ=bΔ,ν=ν,ι=ι,ζ=ζ,x_transform=false)
    result = GenerateSamples!(X, y, R, η=η, nburn=nburn,nsamples=nsamp, V = V, aΔ=aΔ, bΔ=bΔ,ν=ν,ι=ι,ζ=ζ,x_transform=false)
    tock()

    #γ_n = hcat(result.Gammas...)

    γ_n2 = mean(result.Gammas[1:nsamp])
    γ₀ = B₀
    MSE = 0
    for i in 1:190
        MSE = MSE + (γ_n2[i] - γ₀[i])^2
    end
    MSE = MSE * (2/(V*(V-1)))

    return result.Gammas,MSE,result.Xis
end

function output_results(γ::AbstractArray{T},MSE::AbstractFloat,ξ::AbstractArray{U},ξ⁰::DataFrame,simnum,πₛ,μₛ,R,k,jcon::Bool) where {T,U}
    #q = size(γ[1],1)
    #V = convert(Int,(1 + sqrt(1 + 8*q))/2)

    #plot_γ_sim(γ, "Gamma", "simulation$(simnum)_case$(casenum)",jcon)
    show(stdout,"text/plain",DataFrame(MSE=MSE))
    println("")

    nsamp = size(γ[1],1)
    gam = DataFrame(mean=mean(γ))

    γ_sorted = sort.(γ)
    lw = convert(Int64, round(nsamp * 0.05))
    hi = convert(Int64, round(nsamp * 0.95))
    gam[:,"0.05"] = γ_sorted[lw]
    gam[:,"0.95"] = γ_sorted[hi]

    lw = convert(Int64, round(nsamp * 0.025))
    hi = convert(Int64, round(nsamp * 0.975))
    gam[:,"0.025"] = γ_sorted[lw]
    gam[:,"0.975"] = γ_sorted[hi]
    
    gam[:,"pi"] .= πₛ
    gam[:,"mu"] .= μₛ
    gam[:,"R"] .= R

    output = DataFrame(ξ⁰)
    output[:,"Xi posterior"] = ξ
    output[:,"pi"] .= πₛ
    output[:,"mu"] .= μₛ
    output[:,"R"] .= R
    #output[!,"MSE"] = MSE
    mse_df = DataFrame(MSE=MSE)
    mse_df[:,"pi"] .= πₛ
    mse_df[:,"mu"] .= μₛ
    mse_df[:,"R"] .= R

    if jcon
        CSV.write("juliacon/results/sim$(simnum)_pi$(πₛ)-mu$(μₛ)-R$(R)-$(k)microbes_nodes.csv",output)
        CSV.write("juliacon/results/sim$(simnum)_pi$(πₛ)-mu$(μₛ)-R$(R)-$(k)microbes_edges.csv",gam)
        CSV.write("juliacon/results/simu$(simnum)_pi$(πₛ)-mu$(μₛ)-R$(R)-$(k)microbes_MSE.csv",mse_df)
    else
        CSV.write("results/simulation/sim$(simnum)_pi$(πₛ)-mu$(μₛ)-R$(R)-$(k)microbes_nodes.csv",output)
        CSV.write("results/simulation/sim$(simnum)_pi$(πₛ)-mu$(μₛ)-R$(R)-$(k)microbes_edges.csv",gam)
        CSV.write("results/simulation/simu$(simnum)_pi$(πₛ)-mu$(μₛ)-R$(R)-$(k)microbes_MSE.csv",mse_df)
    end
end

#run_case_and_output(3,3,1,0.8,0.1,5,20,false)
#@trace run_case_and_output(10,10,1,1,false)
#@profview run_case_and_output(10,10,1,1,false)
main()
