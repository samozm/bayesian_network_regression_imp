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
        default = 120000
    "--nsamp", "-a"
        help="Number of Gibbs samples to keep (after burn-in)"
        arg_type = Int
        default = 40000
    "--simnum", "-n"
        help="which simulation to run"
        arg_type = Int
        default = 1
    "--samptaxa","-k"
        help = "number of taxa actually used in each sample (relates to filename)"
        arg_type = Int
        default = 15
    "--mean", "-m"
        help = "mean to use to draw edge coefficients (in unrealistic sims) or node coefficients (in realistic sims) - relates to filename"
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
    "--simtype", "-y"
        help = "type of simulation (for realistic sims): additive_phylo, additive_random, interaction_phylo, interaction_random, redundant_phylo, or redundant_random"
        arg_type = String
        default = nothing
        required = false
    "--edgemean", "-e"
        help = "value of the mean of the normal distribution true edge coefficients were drawn from - only used in realistic simulations."
    "--samplesize", "-z"
        help = "sample size (relates to filename)"
        arg_type = Int
        default = 100
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
    type = parsed_CL_args["simtype"]
    edge_μ = parsed_CL_args["edgemean"]
    sampsize = parsed_CL_args["samplesize"]
    Random.seed!(seed)

    run_case_and_output(nburn,nsamp,simnum,μₛ,πₛ,R,k,ν,jcon,type,edge_μ,sampsize)
end

function run_case_and_output(nburn,nsamp,simnum,μₛ,πₛ,R,k,ν,jcon,type="",edge_μ=0.0,sampsize=100)
    loadinfo = Dict("simnum"=>simnum,"pi"=>πₛ,"mu"=>μₛ,"n_microbes"=>k,"out"=>"xis","samplesize"=>sampsize)
    simtypes = Dict(1 => "unrealistic", 2 => "realistic")
    if simtypes[simnum] == "realistic"
        loadinfo["type"] = type
        loadinfo["edge_mu"] = edge_μ
    end
    γ,γ₀,MSE,MSEy,ξ,μ = sim_one_case(nburn,nsamp,loadinfo,jcon,simtypes,simnum,R=R,ν=ν,seed=seed)

    loadinfo["out"] = "xis"
    if jcon
        ξ_in = DataFrame(CSV.File(projectdir("juliacon","data",savename(loadinfo,"csv",digits=1))))
    else
        ξ_in = DataFrame(CSV.File(datadir(joinpath("simulation",simtypes[simnum]),savename(loadinfo,"csv",digits=1))))

        step = Int(floor(nsamp/100))
        if step==0
            step=1
        end
    end

    loadinfo["R"] = R
    loadinfo["nu"] = ν
    output_results(γ[:,:,1],γ₀,MSE,MSEy,mean(ξ,dims=1)[1,:,1],ξ_in,μ,loadinfo,jcon,simtypes,simnum)
end

function sim_one_case(nburn,nsamp,loadinfo,jcon::Bool,simtypes,simnum;η=1.01,ζ=1.0,ι=1.0,R=5,aΔ=1.0,bΔ=1.0,ν=10)
    loadinfo["out"] = "XYs"
    if jcon
        data_in = DataFrame(CSV.File(projectdir("juliacon","data",savename(loadinfo,"csv",digits=1))))

        loadinfo["out"] = "bs"
        b_in = DataFrame(CSV.File(projectdir("juliacon","data",savename(loadinfo,"csv",digits=1))))
        B₀ = convert(Array{Float64,1},b_in[!,:B])
    else
        data_in = DataFrame(CSV.File(datadir(joinpath("simulation",simtypes[simnum]),savename(loadinfo,"csv",digits=1))))

        loadinfo["out"] = "bs"
        b_in = DataFrame(CSV.File(datadir(joinpath("simulation",simtypes[simnum]),savename(loadinfo,"csv",digits=1))))
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

    n = size(y,1)
    y_pred = zeros(n)
    μ_post = mean(result.μ[nburn+1:nburn+nsamp,1,1])
    for i in 1:n
        y_pred[i] = μ_post + sum(γ_n2[1,:,1] .* X[i,:])
    end
    MSEy = sum((y - y_pred).^2)/n

    return result.γ[nburn+1:nburn+nsamp,:,:],γ₀,MSE,MSEy,result.ξ[nburn+1:nburn+nsamp,:,:],result.μ[nburn+1:nburn+nsamp,1,1]
end

function output_results(γ::AbstractArray{T},γ₀::AbstractVector{S},MSE::AbstractFloat,MSEy,ξ::AbstractArray{U},ξ⁰::DataFrame,μ,saveinfo,jcon::Bool,simtypes,simnum,realistic_type="",μₛ=0.0) where {S,T,U}
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
    mse_df[:,"MSEy"] .= MSEy
    mse_df[:,"pi"] .= saveinfo["pi"]
    mse_df[:,"mu"] .= saveinfo["mu"]
    mse_df[:,"R"] .= saveinfo["R"]
    mse_df[:,"n_microbes"] .= saveinfo["n_microbes"]
    mse_df[:,"nu"] .= saveinfo["nu"]

    # this is posterior mu
    μ_sorted = sort(μ)
    mu_df = DataFrame(mean=mean(μ))
    mu_df[:,"0.025"] .= μ_sorted[lw,1,1]
    mu_df[:,"0.975"] .= μ_sorted[hi,1,1]

    mu_df[:,"pi"] .= saveinfo["pi"]
    mu_df[:,"mu"] .= saveinfo["mu"]
    mu_df[:,"R"] .= saveinfo["R"]
    mu_df[:,"n_microbes"] .= saveinfo["n_microbes"]
    mu_df[:,"nu"] .= saveinfo["nu"]

    if simtypes[simnum] == "realistic"
        gam[:,"type"] .= realistic_type
        output[:,"type"] .= realistic_type
        mse_df[:,"type"] .= realistic_type
        mu_df[:,"type"] .= realistic_type
        gam[:,"edge_mu"] .= μₛ
    end

    if jcon
        #TODO better way of adding suffix
        saveinfo["out"] = "nodes"
        CSV.write(projectdir("juliacon","results",savename(saveinfo,"csv",digits=1)),output)
        saveinfo["out"] = "edges"
        CSV.write(projectdir("juliacon","results",savename(saveinfo,"csv",digits=1)),gam)
        saveinfo["out"] = "MSE"
        CSV.write(projectdir("juliacon","results",savename(saveinfo,"csv",digits=1)),mse_df)
    else
        print("out to ")
        println(projectdir("results","simulation",simtypes[simnum]))
        saveinfo["out"] = "nodes"
        CSV.write(projectdir("results","simulation",simtypes[simnum],savename(saveinfo,"csv",digits=1)),output)
        saveinfo["out"] = "edges"
        CSV.write(projectdir("results","simulation",simtypes[simnum],savename(saveinfo,"csv",digits=1)),gam)
        saveinfo["out"] = "MSE"
        CSV.write(projectdir("results","simulation",simtypes[simnum],savename(saveinfo,"csv",digits=1)),mse_df)
        saveinfo["out"] = "mu"
        CSV.write(projectdir("results","simulation",simtypes[simnum],savename(saveinfo,"csv",digits=1)),mu_df)
    end
end

#@profview run_case_and_output(100,100,1,0.8,0.3,9,22,10,false,"",0.0,500)
main()
