using LinearAlgebra

using CSV,ArgParse

using Random, DataFrames, StatsBase, InvertedIndices, ProgressMeter, Distributions
using StaticArrays,TypedTables
using BayesianNetworkRegression,DrWatson,MCMCDiagnosticTools,JLD2
#include("../BayesianNetworkRegression.jl/src/gelmandiag.jl")

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
    "--timeout", "-t"
        help = "flag indicating time should be saved"
        action = :store_true
    end
    return parse_args(args)
end

function main()
    parsed_CL_args = parse_CL_args()
    nburn = parsed_CL_args["nburn"]
    nsamp = parsed_CL_args["nsamp"]
    simnum = parsed_CL_args["simnum"]
    μₛ = parsed_CL_args["mean"]
    πₛ = parsed_CL_args["pi"]
    R = parsed_CL_args["r"]
    k = parsed_CL_args["samptaxa"]
    ν = parsed_CL_args["nu"]
    seed = parsed_CL_args["seed"]
    typ = parsed_CL_args["simtype"]
    edge_μ = parsed_CL_args["edgemean"]
    sampsize = parsed_CL_args["samplesize"]
    tmot = parsed_CL_args["timeout"]
    Random.seed!(seed)

    run_case_and_output(nburn,nsamp,simnum,μₛ,πₛ,R,k,ν,typ,edge_μ,sampsize,seed,tmot)
end

function run_case_and_output(nburn,nsamp,simnum,μₛ,πₛ,R,k,ν,typ="",edge_μ=0.0,sampsize=100,seed=nothing,tmot=false)
    loadinfo = Dict("simnum"=>simnum,"pi"=>πₛ,"mu"=>μₛ,"n_microbes"=>k,"out"=>"xis","samplesize"=>sampsize)
    simtypes = Dict(1 => "unrealistic", 2 => "realistic")
    if simtypes[simnum] == "realistic"
        loadinfo["type"] = typ
        loadinfo["edge_mu"] = edge_μ

        if μₛ==0.8 && πₛ==0.3 && k==8 && (typ=="additive_random" || typ=="additive_phylo") && sampsize==1000
            nburn = 200000
        elseif μₛ==1.6 && πₛ==0.8 && k==22 && (typ=="redundant_random" || typ=="redundant_phylo") && sampsize==500
            nburn = 160000
        elseif μₛ==0.8 && πₛ==0.3 && k==22 && typ=="redundant_phylo" && sampsize==500
            nburn = 200000
        elseif μₛ==0.8 && πₛ==0.3 && k==22 && typ=="redundant_phylo" && sampsize==1000
            nburn = 300000
        end
    end
    
    sim_one_case(nburn,nsamp,loadinfo,simtypes,simnum,seed=seed,η=1.01,ζ=1.0,ι=1.0,R=R,aΔ=1.0,bΔ=1.0,ν=ν,tmot=tmot)
end

function sim_one_case(nburn,nsamp,loadinfo,simtypes,simnum;seed=nothing,η=1.01,ζ=1.0,ι=1.0,R=5,aΔ=1.0,bΔ=1.0,ν=10,tmot=false)
    loadinfo["out"] = "XYs"
    #data_in = DataFrame(CSV.File(projectdir(simtypes[simnum],"data",savename(loadinfo,"csv",digits=1))))
    data_in = DataFrame(CSV.File(string(simtypes[simnum],"/data/",savename(loadinfo,"csv",digits=1))))

        

    X = Matrix(data_in[:,names(data_in,Not("y"))])
    y = SVector{size(X,1)}(data_in[:,:y])

    q = size(X,2)
    V = convert(Int,(1 + sqrt(1 + 8*q))/2)
    
    try 
       @show size(X), size(y), R, η, nburn,nsamp, V, aΔ, bΔ,ν,ι,ζ,num_chains, seed
    catch
    end 

    tm=@elapsed result = Fit!(X, y, R, η=η, nburn=nburn,nsamples=nsamp, aΔ=aΔ, bΔ=bΔ,ν=ν,ι=ι,ζ=ζ,x_transform=false,num_chains=num_chains,seed=seed,in_seq=true)
    
    loadinfo["out"] = "bs"
    b_in = DataFrame(CSV.File(string(simtypes[simnum],"/data/",savename(loadinfo,"csv",digits=1))))
    B₀ = convert(Array{Float64,1},b_in[!,:B])

    loadinfo["out"] = "xis"

    ξ_in = DataFrame(CSV.File(string(simtypes[simnum],"/data/",savename(loadinfo,"csv",digits=1))))

    loadinfo["R"] = R
    loadinfo["nu"] = ν

    γ_n2 = mean(result.state.γ[:,:,:],dims=1)
    γ₀ = B₀
    MSE = 0

    for i in 1:size(γ_n2,2)
        MSE = MSE + (γ_n2[i] - γ₀[i])^2
    end
    MSE = MSE * (2/(V*(V-1)))

    n = size(y,1)
    y_pred = zeros(n)
    μ_post = mean(result.state.μ[:,1,1])
    for i in 1:n
        y_pred[i] = μ_post + sum(γ_n2[1,:,1] .* X[i,:])
    end
    MSEy = sum((y - y_pred).^2)/n
    output_results(result.state.γ, γ₀, MSE, MSEy,mean(result.state.ξ,dims=1)[1,:,1],ξ_in,result.state.μ[:,1,1],result.psrf,loadinfo,simtypes,simnum)
    time_df = DataFrame(time=tm)
    time_df[:,"pi"] .= loadinfo["pi"]
    time_df[:,"mu"] .= loadinfo["mu"]
    time_df[:,"R"] .= loadinfo["R"]
    time_df[:,"n_microbes"] .= loadinfo["n_microbes"]
    time_df[:,"nu"] .= loadinfo["nu"]
    loadinfo["out"] = "time"
    CSV.write(string(savename(loadinfo,"csv",digits=1)),time_df)
end


function output_results(γ::AbstractArray{T},γ₀::AbstractVector{S},MSE::AbstractFloat,MSEy,ξ::AbstractArray{U},
    ξ⁰::DataFrame,μ,psrf,saveinfo,simtypes,simnum,realistic_type="",μₛ=0.0) where {S,T,U}
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

    psrf_df = DataFrame(psrf)
    
    @show psrf_df

    psrf_df[:,"mean_xi"] .= mean(psrf.ξ[1:V])
    psrf_df[:,"max_xi"] .= max(psrf.ξ[1:V]...)
    psrf_df[:,"mean_gamma"] .= mean(psrf.γ[1:q])
    psrf_df[:,"max_gamma"] .= max(psrf.γ[1:q]...)

    type = simtypes[simnum]
    saveinfo["out"] = "nodes"
    CSV.write(string(type,"-results/",savename(saveinfo,"csv",digits=1)),output)
    saveinfo["out"] = "edges"
    CSV.write(string(type,"-results/",savename(saveinfo,"csv",digits=1)),gam)
    saveinfo["out"] = "MSE"
    CSV.write(string(type,"-results/",savename(saveinfo,"csv",digits=1)),mse_df)
    saveinfo["out"] = "mu"
    CSV.write(string(type,"-results/",savename(saveinfo,"csv",digits=1)),mu_df)
    saveinfo["out"] = "psrf"
    CSV.write(string(type,"-results/",savename(saveinfo,"csv",digits=1)),psrf_df)
    
end


main()
