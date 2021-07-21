using LinearAlgebra
using Base: Bool, Float16, Int16
using CSV,ArgParse,TickTock
using ProfileView, Traceur
#using BayesianNetworkRegression

using DataFrames: Vector
using Core: Typeof
using Base: Float64
using Random, DataFrames, StatsBase, InvertedIndices, ProgressMeter, Distributions
using StaticArrays,TypedTables
using BayesianNetworkRegression
using Random, DataFrames, LinearAlgebra, StatsBase
using Distributions
using BenchmarkTools,ProfileView
using DrWatson

function parse_CL_args()
    args = ArgParseSettings()
    @add_arg_table! args begin
    "--nburn", "-b"
        help="Number of burn-in Gibbs samples to take"
        arg_type = Int
        default = 30000
    "--nsamp", "-a"
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
        help = "R value, the dimension of the latent space (relates to filename)"
        arg_type = Int
        default = 5
    "--juliacon", "-j"
        help = "flag indicating output should go to juliacon folder"
        action = :store_true
    "--seed", "-s"
        help = "random seed to use for simulations"
        arg_type = Int
        default = 734
    "--nu", "-u"
        help = "nu value (hyperparameter)"
        arg_type = Int
        default = 10
    end
    return parse_args(args)
end

function main()
    @quickactivate
    parsed_CL_args = parse_CL_args()
    nburn = parsed_CL_args["nburn"]
    nsamp = parsed_CL_args["nsamp"]
    simnum = parsed_CL_args["simnum"]
    jcon = parsed_CL_args["juliacon"]
    μₛ = parsed_CL_args["mean"]
    πₛ = parsed_CL_args["pi"]
    R = parsed_CL_args["r"]
    k = parsed_CL_args["samptaxa"]
    ν = parsed_CL_args["nu"]
    seed = parsed_CL_args["seed"]
    Random.seed!(seed)

    run_case_and_output(nburn,nsamp,simnum,μₛ,πₛ,R,k,ν,jcon)
end

function run_case_and_output(nburn,nsamp,simnum,μₛ,πₛ,R,k,ν,jcon)
    loadinfo = Dict("simnum"=>simnum,"pi"=>πₛ,"mu"=>μₛ,"n_microbes"=>k,"out"=>"xis")
    γ,γ₀,MSE,ξ = sim_one_case(nburn,nsamp,loadinfo,jcon,R=R,ν=ν)

    loadinfo["out"] = "xis"
    if jcon
        ξ_in = DataFrame(CSV.File(projectdir("juliacon","data",savename(loadinfo,"csv",digits=1))))
    else
        ξ_in = DataFrame(CSV.File(datadir("simulation",savename(loadinfo,"csv",digits=1))))

        step = Int(floor(nsamp/100))
        if step==0
            step=1
        end
    end

    loadinfo["R"] = R
    loadinfo["nu"] = ν
    output_results(γ[:,:,1],γ₀,MSE,mean(ξ,dims=1)[1,:,1],ξ_in,loadinfo,jcon)
end

function sim_one_case(nburn,nsamp,loadinfo,jcon::Bool;η=1.01,ζ=1.0,ι=1.0,R=5,aΔ=1.0,bΔ=1.0,ν=10)
    loadinfo["out"] = "XYs"
    if jcon
        data_in = DataFrame(CSV.File(projectdir("juliacon","data",savename(loadinfo,"csv",digits=1))))

        loadinfo["out"] = "bs"
        b_in = DataFrame(CSV.File(projectdir("juliacon","data",savename(loadinfo,"csv",digits=1))))
        B₀ = convert(Array{Float64,1},b_in[!,:B])
    else
        data_in = DataFrame(CSV.File(datadir("simulation",savename(loadinfo,"csv",digits=1))))

        loadinfo["out"] = "bs"
        b_in = DataFrame(CSV.File(datadir("simulation",savename(loadinfo,"csv",digits=1))))
        B₀ = convert(Array{Float64,1},b_in[!,:B])
    end
    #X = convert(Matrix,data_in[:,names(data_in,Not("y"))])
    X = Matrix(data_in[:,names(data_in,Not("y"))])
    y = SVector{size(X,1)}(data_in[:,:y])

    q = size(X,2)
    V = convert(Int,(1 + sqrt(1 + 8*q))/2)
    
    tick()
    #τ², u, ξ, γ, D, θ, Δ, M, μ, Λ, πᵥ = BayesNet(X, y, R, η=η, nburn=nburn,nsamples=nsamp, V_in = V, aΔ=aΔ, bΔ=bΔ,ν=ν,ι=ι,ζ=ζ,x_transform=false)
    result = GenerateSamples!(X, y, R, η=η, nburn=nburn,nsamples=nsamp, V = V, aΔ=aΔ, bΔ=bΔ,ν=ν,ι=ι,ζ=ζ,x_transform=false)
    tock()

    #γ_n = hcat(result.Gammas...)

    γ_n2 = mean(result.γ[nburn+1:nburn+nsamp,:,:],dims=1)
    γ₀ = B₀
    MSE = 0
    for i in 1:190
        MSE = MSE + (γ_n2[i] - γ₀[i])^2
    end
    MSE = MSE * (2/(V*(V-1)))

    return result.γ[nburn+1:nburn+nsamp,:,:],γ₀,MSE,result.ξ[nburn+1:nburn+nsamp,:,:]
end

function output_results(γ::AbstractArray{T},γ₀::AbstractVector{S},MSE::AbstractFloat,ξ::AbstractArray{U},ξ⁰::DataFrame,saveinfo,jcon::Bool) where {S,T,U}
    q = size(γ,2)
    V = convert(Int,(1 + sqrt(1 + 8*q))/2)

    show(stdout,"text/plain",DataFrame(MSE=MSE))
    println("")

    nsamp = size(γ,1)
    gam = DataFrame(mean=mean(γ,dims=1)[1,:])

    γ_sorted = sort(γ,dims=1)
    lw = convert(Int64, round(nsamp * 0.05))
    hi = convert(Int64, round(nsamp * 0.95))
    if lw == 0
        lw = 1
    end
    gam[:,"0.05"] = γ_sorted[lw,:,1]
    gam[:,"0.95"] = γ_sorted[hi,:,1]

    lw = convert(Int64, round(nsamp * 0.025))
    hi = convert(Int64, round(nsamp * 0.975))
    if lw == 0
        lw = 1
    end
    gam[:,"0.025"] = γ_sorted[lw,:,1]
    gam[:,"0.975"] = γ_sorted[hi,:,1]

    gam[:,"pi"] .= saveinfo["pi"]
    gam[:,"mu"] .= saveinfo["mu"]
    gam[:,"R"] .= saveinfo["R"]
    gam[:,"true_B"] = γ₀
    gam[:,"n_microbes"] .= saveinfo["n_microbes"]

    gam[:,"nu"] .= saveinfo["nu"]

    gam[:,"y_microbe"] .= 0
    gam[:,"x_microbe"] .= 0

    l = 1
    for i in 1:V
        for j in i+1:V
            gam[l,"y_microbe"] = i
            gam[l,"x_microbe"] = j
            l += 1
        end
    end

    output = DataFrame(ξ⁰)
    output[:,"Xi posterior"] = ξ
    output[:,"pi"] .= saveinfo["pi"]
    output[:,"mu"] .= saveinfo["mu"]
    output[:,"R"] .= saveinfo["R"]
    output[:,"n_microbes"] .= saveinfo["n_microbes"]
    output[:,"nu"] .= saveinfo["nu"]

    mse_df = DataFrame(MSE=MSE)
    mse_df[:,"pi"] .= saveinfo["pi"]
    mse_df[:,"mu"] .= saveinfo["mu"]
    mse_df[:,"R"] .= saveinfo["R"]
    mse_df[:,"n_microbes"] .= saveinfo["n_microbes"]
    mse_df[:,"nu"] .= saveinfo["nu"]

    if jcon
        #TODO better way of adding suffix
        saveinfo["out"] = "nodes"
        CSV.write(projectdir("juliacon","results",savename(saveinfo,"csv",digits=1)),output)
        saveinfo["out"] = "edges"
        CSV.write(projectdir("juliacon","results",savename(saveinfo,"csv",digits=1)),gam)
        saveinfo["out"] = "MSE"
        CSV.write(projectdir("juliacon","results",savename(saveinfo,"csv",digits=1)),mse_df)
    else
        saveinfo["out"] = "nodes"
        CSV.write(projectdir("results","simulation",savename(saveinfo,"csv",digits=1)),output)
        saveinfo["out"] = "edges"
        CSV.write(projectdir("results","simulation",savename(saveinfo,"csv",digits=1)),gam)
        saveinfo["out"] = "MSE"
        CSV.write(projectdir("results","simulation",savename(saveinfo,"csv",digits=1)),mse_df)
    end
end

#cd("bayesian_network_regression_imp")
#@quickactivate
#run_case_and_output(3,3,1,0.8,0.3,5,15,false)
#@trace run_case_and_output(10,10,1,1,false)
#@profview run_case_and_output(3,3,1,0.8,0.1,5,20,false)
main()
