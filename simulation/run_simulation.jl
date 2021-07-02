using LinearAlgebra: Matrix
using Base: Bool, Float16, Int16
using CSV,ArgParse,TickTock
using ProfileView, Traceur
using BayesianNetworkRegression
using Random, DataFrames, LinearAlgebra, StatsBase
using Distributions
using BenchmarkTools
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
    "--casenum", "-c"
        help="which case of the simulation to run"
        arg_type = Int
        default = 1
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
    casenum = parsed_CL_args["casenum"]
    jcon = parsed_CL_args["juliacon"]

    run_case_and_output(nburn,nsamp,simnum,casenum,jcon)
end

function run_case_and_output(nburn,nsamp,simnum,casenum,jcon)
    γ,MSE,ξ = sim_one_case(simnum,casenum,nburn,nsamp,jcon)

    if jcon
        ξ_in = DataFrame(CSV.File("juliacon/data/simulation$(simnum)_case$(casenum)_xis.csv"))
    else
        ξ_in = DataFrame(CSV.File("data/simulation/simulation$(simnum)_case$(casenum)_xis.csv"))

        step = Int(floor(nsamp/100))
        if step==0
            step=1
        end
        ξ_prog = map(i->mean(ξ[1:i]),1:step:nsamp)
        #println("tst")
        #show(stdout,"text/plain",tst)
        #println("")
        savefig(plot(transpose(hcat(ξ_prog...)),legend=false),"plots/simulation/xi_converge_simulation$(simnum)_case$(casenum)")

        γ_prog = zeros(190,length(1:step:nsamp))
        for j in 1:190
            γ_prog[j,:] = map(i->mean(γ[j,1:i]),1:step:nsamp)
        end
        savefig(plot(transpose(γ_prog),legend=false),"plots/simulation/gamma_coverage_simulation$(simnum)_case$(casenum)")

    end

    output_results(γ[:,1:nsamp],MSE,mean(ξ[1:nsamp]),ξ_in,simnum,casenum,jcon)
end

function sim_one_case(simnum::Int64,casenum::Int64,nburn::Int64,nsamp::Int64,jcon::Bool,η::Float64=1.01,ζ::Float64=1.0,ι::Float64=1.0,R::Int64=5,aΔ::Float64=1.0,bΔ::Float64=1.0,ν::Int64=10)
    if jcon
        data_in = DataFrame(CSV.File("juliacon/data/simulation$(simnum)_case$(casenum).csv"))

        b_in = DataFrame(CSV.File("juliacon/data/simulation$(simnum)_case$(casenum)_bs.csv"))
        B₀ = convert(Array{Float64,1},b_in[!,:B])
    else
        data_in = DataFrame(CSV.File("data/simulation/simulation$(simnum)_case$(casenum).csv"))

        b_in = DataFrame(CSV.File("data/simulation/simulation$(simnum)_case$(casenum)_bs.csv"))
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

    low = zeros(190)
    high = zeros(190)
    lw = convert(Int64, round(nsamp * 0.05))
    hi = convert(Int64, round(nsamp * 0.95))

    γ_n = hcat(result.Gammas...)

    γ_n2 = mean(result.Gammas[1:nsamp])
    γ₀ = B₀
    MSE = 0
    for i in 1:190
        MSE = MSE + (γ_n2[i] - γ₀[i])^2
    end
    MSE = MSE * (2/(V*(V-1)))

    return γ_n,MSE,result.Xis
end

function output_results(γ::Array,MSE::Float64,ξ::Array,ξ⁰::DataFrame,simnum::Int64,casenum::Int64,jcon::Bool)
    q = size(γ,1)
    V = convert(Int,(1 + sqrt(1 + 8*q))/2)

    plot_γ_sim(γ, "Gamma", "simulation$(simnum)_case$(casenum)",jcon)
    #ξ⁰[!,"Xi Probability"] = ξ
    #show(stdout,"text/plain",ξ⁰)
    #println("")
    show(stdout,"text/plain",DataFrame(MSE=MSE))
    println("")
    #println("Gamma")
    #gam = DataFrame(create_upper_tri(vec(median.(γ)),V),:auto)
    gam = DataFrame(create_upper_tri(vec(mean.(γ)),V),:auto)
    #show(stdout,"text/plain",gam[!,names(gam,Not("x1"))])
    #println("")
    output = DataFrame(ξ⁰)
    output[:,"Xi posterior"] = ξ
    #output[!,"MSE"] = MSE
    if jcon
        CSV.write("juliacon/results/simulation$(simnum)_case$(casenum).csv",output)
        CSV.write("juliacon/results/simulation$(simnum)_case$(casenum)_gammas.csv",gam)
    else
        CSV.write("results/simulation/simulation$(simnum)_case$(casenum).csv",output)
        CSV.write("results/simulation/simulation$(simnum)_case$(casenum)_gammas.csv",gam)
        CSV.write("results/simulation/simulation$(simnum)_case$(casenum)_MSE.csv",DataFrame(MSE=MSE))
    end
end

run_case_and_output(3,3,1,1,false)
#@trace run_case_and_output(10,10,1,1,false)
#main()
