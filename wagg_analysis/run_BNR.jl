using DrWatson
@quickactivate

using LinearAlgebra
using CSV,ArgParse
using Random, DataFrames, StatsBase, InvertedIndices, ProgressMeter, Distributions
using StaticArrays,TypedTables
using BayesianNetworkRegression,DrWatson,MCMCDiagnosticTools,JLD2,Distributed,DrWatson

MIN_GEN = 10000 #10000 #20000
MAX_GEN = 200000 #200000 #800000

NUM_CHAINS = 3

NORMED = false

addprocs(NUM_CHAINS,exeflags=["--optimize=0","--math-mode=ieee","--check-bounds=yes"])

@everywhere begin
    #using DrWatson
    #@quickactivate
    #using BayesianNetworkRegression,CSV,DataFrames,StaticArrays
    include("../../BayesianNetworkRegression.jl/src/BayesianNetworkRegression.jl")
    using CSV,DataFrames,StaticArrays
    using TypedTables,Random,LinearAlgebra,Distributions,DrWatson
end

function main()
    seed = 102722
    sampsize = 100
    cutoff = "40"
    tmot = true
    Random.seed!(seed)

    ### SESH 1
    #Rs    = [   5,   5,   5,    6,   6,   6]
    #nu    = [   9,   7,   8,   10,   8,   9]
    #Rs    = [ 5,  6,  7,  7]
    #nu    = [12, 12, 10, 11]

    #Rs    = [ 9, 10, 11, 12]
    #nu    = [11, 13, 14, 15]

    #Rs    = [15, 15, 15, 15]
    #nu    = [17, 18, 19, 20]

    # Done (50)
    #Rs    = [13, 13, 15, 15]
    #nu    = [16, 17, 18, 19]

    #Rs    = [5, 5, 7,  7]
    #nu    = [7, 9, 9, 11]

    # Done
    #Rs    = [6,  6,  8,  8]
    #nu    = [8, 10, 10, 12]

    Rs    = [20, 20,  8,  6]
    nu    = [24, 26, 16, 12]

    Rs    = [30, 20, 40, 50]
    nu    = [36, 36, 46, 56]

    Rs    = [6,   6,  6]
    nu    = [18, 20, 24]

    #Rs    = [4,   6, 24, 24]
    #nu    = [14, 16, 28, 30]

    # Done
    #Rs    = [2, 2, 10, 10]
    #nu    = [4, 6, 12, 14]

    ns = ["200"] #["100","200"]

    #Rs    = [ 5,  7,  8,  9,  8]
    #nu    = [12, 12, 12, 12, 11]
    nburn = 5000
    nsamp = 5000

    idx = 1
    n_cases = size(Rs)[1]
    for n in ns
        for i in 1:n_cases
            println(stderr,"-------------------------")
            println(stderr,"case ", idx, " of ", n_cases*size(ns,1))
            println(stderr,"R=",Rs[i]," nu=",nu[i]," n=",n)
            idx = idx+1
            run_case_and_output(nburn,nsamp,cutoff,Rs[i],nu[i],n,sampsize,seed,tmot)
            GC.gc()
        end
    end
end

function run_case_and_output(nburn,nsamp,cutoff,R,ν,nsamples,sampsize=100,seed=nothing,tmot=false)
    loadinfo = Dict("cutoff"=>cutoff)#,"samplesize"=>sampsize)
    
    sim_one_case(nburn,nsamp,loadinfo,nsamples,seed=seed,η=1.01,ζ=1.0,ι=1.0,R=R,aΔ=1.0,bΔ=1.0,ν=ν,tmot=tmot)
end

function sim_one_case(nburn,nsamp,loadinfo,nsamples;seed=nothing,η=1.01,ζ=1.0,ι=1.0,R=5,aΔ=1.0,bΔ=1.0,ν=10,tmot=false)
    loadinfo["out"] = "XY"
    loadinfo["n"] = nsamples
    data_in = DataFrame(CSV.File(string("bayesian_network_regression_imp/data/wagg/",savename(loadinfo,"csv",digits=1))))


    X = Matrix(data_in[:,names(data_in,Not(["Nleach","Pleach","Decom"]))])
    y = SVector{size(X,1)}(data_in[:,:Pleach])
    x_names = names(data_in,Not(["Nleach","Pleach","Decom"]))

    q = size(X,2)
    V = convert(Int,round((-1 + sqrt(1 + 8*q))/2))
    n = size(X,1)

    if NORMED
        onefill = max(X[X .!= 1 .&& X .!=0]...) #mean(X[X .!= 1 .&& X .!=0])

        for i in 1:n
            for j in 1:q
                if X[i,j] == 1
                    X[i,j] = onefill
                end
            end
        end
        #println(stderr,"onefill=",onefill)
        #println(stderr,"size(X)=",size(X))
        #println(stderr,"X[1:10,1:10]")
        #println(stderr, X[1:10,1:10])
    end

    @everywhere begin
        num_chains = $(NUM_CHAINS)
        R = $(R)
        V = $(V)
        mingen = $(MIN_GEN)
        maxgen = $(MAX_GEN)
        nburn = $(nburn)
        nsamp = $(nsamp)
        total = nburn+nsamp
        q = floor(Int,V*(V-1)/2)
        seed = $(seed)

        loadinfo = $(loadinfo)


        Random.seed!(seed)

        X = $(X)
        y = $(y)
    end
    num_chains=3
    purge_burn=1000
    #purge_burn=50

    num_chains = NUM_CHAINS
    psrf_cutoff = 1.01

    println("start fit")
    #tm=@elapsed result = Fit!(X, y, R, η=η, nburn=nburn,nsamples=nsamp,aΔ=aΔ, 
    #                            bΔ=bΔ,ν=ν,ι=ι,ζ=ζ,x_transform=false,suppress_timer=false,
    #                            num_chains=num_chains,seed=seed,
    ##                            purge_burn=purge_burn)
    tm=@elapsed result = BayesianNetworkRegression.Fit!(X, y, R, η=η, V=V, mingen=mingen, 
                                    maxgen=maxgen, aΔ=aΔ, 
                                    bΔ=bΔ,ν=ν,ι=ι,ζ=ζ,x_transform=false,suppress_timer=false,
                                    num_chains=num_chains,seed=seed,
                                    purge_burn=purge_burn,psrf_cutoff=psrf_cutoff)

    # This must be done before we re-assign nburn
    loadinfo["R"] = string(R)
    loadinfo["nu"] = string(ν)
    #loadinfo["nburn"] = string(nburn)
    #loadinfo["nsamples"] = string(nsamp)

    nburn = convert(Int64, round(MIN_GEN / 2))

    if purge_burn < nburn
        if nburn % purge_burn != 0 
            nburn = purge_burn - (nburn % purge_burn)
        else
            nburn = purge_burn
        end
    end

    total = size(result.state.γ,1)

    if NORMED
        loadinfo["normed"] = "yes"
    end

    
    output_results(result.state.γ[nburn+1:total,:,:],mean(result.state.ξ[nburn+1:total,:,:],dims=1)[1,:,1],
                   result.state.μ[nburn+1:total,1,1],result.rhatξ,result.rhatγ,loadinfo,x_names)
    time_df = DataFrame(time=tm)
    time_df[:,"cutoff"] .= loadinfo["cutoff"]
    loadinfo["out"] = "time"
    time_df[:,"tot_samples"] .= result.sampled * 2
    println("")
    println("-----------------------------------------------")
    println("")
    CSV.write(string("results/wagg/",savename(loadinfo,"csv",digits=1)),time_df)
end


function output_results(γ::AbstractArray{T},ξ::AbstractArray{U},μ,rhatξ,rhatγ,saveinfo,names) where {T,U}
    q = size(γ,2)
    V = convert(Int,(-1 + sqrt(1 + 8*q))/2)


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

    gam[:,"cutoff"] .= saveinfo["cutoff"]

    gam[:,"microbes"] .= names

    #l = 1
    #for i in 1:V
    #    for j in i+1:V
    #        gam[l,"y_microbe"] = names[i]
    #        gam[l,"x_microbe"] = names[j]
    #        l += 1
    #    end
    #end
    mics_mat = Matrix{Union{Missing, String}}(missing,V,V)
    i = 1
    for k = 1:V
        for l = k:V
            mics_mat[l,k] = split(names[i,1],"x")[1]
            i += 1
        end
    end

    mics = mics_mat[V,:]
    mics[size(mics,1)] = split(names[q,1],"x")[2]

    output = DataFrame(Xi_posterior=ξ)
    output[:,"cutoff"] .= saveinfo["cutoff"]
    output[:,"microbe"] = mics

    # this is posterior mu
    μ_sorted = sort(μ)
    mu_df = DataFrame(mean=mean(μ))
    mu_df[:,"0.025"] .= μ_sorted[lw,1,1]
    mu_df[:,"0.975"] .= μ_sorted[hi,1,1]

    mu_df[:,"cutoff"] .= saveinfo["cutoff"]
    

    psrf_df = DataFrame(rhatγ=rhatγ.γ)
    psrf_df[:,"rhatξ"] = vcat(rhatξ.ξ,zeros(q-V))
    
    psrf_df[:,"mean_xi"] .= mean(rhatξ.ξ[1:V,])
    psrf_df[:,"max_xi"] .= max(rhatξ.ξ[1:V,]...)
    psrf_df[:,"mean_gamma"] .= mean(rhatγ.γ[1:q,])
    psrf_df[:,"max_gamma"] .= max(rhatγ.γ[1:q,]...)

    @show psrf_df[1,"max_xi"];@show psrf_df[1,"max_gamma"]
    
    saveinfo["out"] = "nodes"
    CSV.write(string("results/wagg/",savename(saveinfo,"csv",digits=1)),output)

    saveinfo["out"] = "edges"
    CSV.write(string("results/wagg/",savename(saveinfo,"csv",digits=1)),gam)

    saveinfo["out"] = "mu"
    CSV.write(string("results/wagg/",savename(saveinfo,"csv",digits=1)),mu_df)
    
    saveinfo["out"] = "psrf"
    CSV.write(string("results/wagg/",savename(saveinfo,"csv",digits=1)),psrf_df)
    
end


main()
