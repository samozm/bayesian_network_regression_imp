using LinearAlgebra

using CSV,ArgParse

using Random, DataFrames, StatsBase, InvertedIndices, ProgressMeter, Distributions
using StaticArrays,TypedTables
using BayesianNetworkRegression,DrWatson,MCMCDiagnosticTools,JLD2,Distributed

addprocs(3,exeflags=["--optimize=0","--math-mode=ieee","--check-bounds=yes"])

@everywhere begin
    using BayesianNetworkRegression,CSV,DataFrames,StaticArrays
    using TypedTables,Random,LinearAlgebra,Distributions,DrWatson
end

function main()
    #nburn = 120000
    #nsamp = 40000
    simnum = 2
    μₛs = [0.8,1.6]
    πₛs = [0.3,0.8]
    R = 7
    edge_μ = 0.4
    ks = [8,22]
    ks = [8,22]
    ν = 10
    seed = 2358
    sampsizes = [500,1000]
    simtypes = ["additive_random","additive_phylo","interaction_random","interaction_phylo","redundant_random","redundant_phylo"]
    tmot = true
    Random.seed!(seed)

    #nburn = 60000
    #nsamp = 40000
    nburn = 30000
    nsamp = 20000

    for μₛ in μₛs
        for πₛ in πₛs
            for k in ks
                for sampsize in sampsizes
                    for typ in simtypes
                        #if (μₛ==1.6 && πₛ==0.3 && k==8 && sampsize==500 && typ=="additive_random")
                        #    nburn = 120000 #prev 60000
                        #elseif (μₛ==1.6 && πₛ==0.3 && k==22 && sampsize==500 && typ=="redundant_phylo")
                        #    nburn = 100000 #prev 60000
                        #elseif (μₛ==1.6 && πₛ==0.3 && k==22 && sampsize==1000 && typ=="redundant_random")
                        #    nburn = 60000
                        if (μₛ==0.8 && πₛ==0.8 && k==22 && sampsize==1000 && typ=="additive_random")
                            #nburn = 60000
                        else continue
                        end
                        run_case_and_output(nburn,nsamp,simnum,μₛ,πₛ,R,k,ν,typ,edge_μ,sampsize,seed,tmot)
                        #nburn = 60000
                        nburn = 30000
                        GC.gc()
                    end
                end
            end
        end
    end
end

function run_case_and_output(nburn,nsamp,simnum,μₛ,πₛ,R,k,ν,typ="",edge_μ=0.0,sampsize=100,seed=nothing,tmot=false)
    loadinfo = Dict("simnum"=>simnum,"pi"=>πₛ,"mu"=>μₛ,"n_microbes"=>k,"out"=>"xis","samplesize"=>sampsize)
    simtypes = Dict(1 => "unrealistic", 2 => "realistic")
    loadinfo["type"] = typ
    loadinfo["edge_mu"] = edge_μ
    
    @show nburn,nsamp
    @show μₛ; @show πₛ; @show R; @show k; @show ν; @show sampsize; @show typ
    sim_one_case(nburn,nsamp,loadinfo,simtypes,simnum,seed=seed,η=1.01,ζ=1.0,ι=1.0,R=R,aΔ=1.0,bΔ=1.0,ν=ν,tmot=tmot)
end

function sim_one_case(nburn,nsamp,loadinfo,simtypes,simnum;seed=nothing,η=1.01,ζ=1.0,ι=1.0,R=5,aΔ=1.0,bΔ=1.0,ν=10,tmot=false)
    loadinfo["out"] = "XYs"
    data_in = DataFrame(CSV.File(string("bayesian_network_regression_imp/data/simulation/",simtypes[simnum],"/",savename(loadinfo,"csv",digits=1))))

        

    X = Matrix(data_in[:,names(data_in,Not("y"))])
    y = SVector{size(X,1)}(data_in[:,:y])

    q = size(X,2)
    V = convert(Int,(1 + sqrt(1 + 8*q))/2)


    if simnum == 1
        println("wrong simnum")
    else
        @everywhere begin
        
            R = $(R)
            V = $(V)
            nburn = $(nburn)
            nsamp = $(nsamp)
            total = nburn+nsamp
            q = floor(Int,V*(V-1)/2)
            seed = $(seed)
            simnum=1

            loadinfo = $(loadinfo)
            simtypes = Dict(1 => "unrealistic", 2 => "realistic")


            Random.seed!(seed)

            X = $(X)
            y = $(y)
        end
        num_chains=3
        purge_burn=10000
        tm=@elapsed result = Fit!(X, y, R, η=η, nburn=nburn,nsamples=nsamp, aΔ=aΔ, 
                                    bΔ=bΔ,ν=ν,ι=ι,ζ=ζ,x_transform=false,suppress_timer=false,
                                    num_chains=num_chains,seed=seed,
                                    purge_burn=purge_burn)

        nburn = purge_burn
        

        loadinfo["out"] = "bs"
        b_in = DataFrame(CSV.File(string("bayesian_network_regression_imp/data/simulation/",simtypes[simnum],"/",savename(loadinfo,"csv",digits=1))))
        B₀ = convert(Array{Float64,1},b_in[!,:B])

        loadinfo["out"] = "xis"
        ξ_in = DataFrame(CSV.File(string("bayesian_network_regression_imp/data/simulation/",simtypes[simnum],"/",savename(loadinfo,"csv",digits=1))))

        total = nburn + nsamp
        loadinfo["R"] = R
        loadinfo["nu"] = ν

        γ_n2 = mean(result.state.γ[nburn+1:total,:,:],dims=1)
        γ₀ = B₀
        MSE = 0

        for i in 1:size(γ_n2,2)
            MSE = MSE + (γ_n2[i] - γ₀[i])^2
        end
        MSE = MSE * (2/(V*(V-1)))

        n = size(y,1)
        y_pred = zeros(n)
        μ_post = mean(result.state.μ[nburn+1:total,1,1])
        for i in 1:n
            y_pred[i] = μ_post + sum(γ_n2[1,:,1] .* X[i,:])
        end
        MSEy = sum((y - y_pred).^2)/n
        output_results(result.state.γ[nburn+1:total,:,:], γ₀, MSE, MSEy,mean(result.state.ξ[nburn+1:total,:,:],dims=1)[1,:,1],ξ_in,result.state.μ[nburn+1:total,1,1],result.rhat,loadinfo,simtypes,simnum)
        time_df = DataFrame(time=tm)
        time_df[:,"pi"] .= loadinfo["pi"]
        time_df[:,"mu"] .= loadinfo["mu"]
        time_df[:,"R"] .= loadinfo["R"]
        time_df[:,"n_microbes"] .= loadinfo["n_microbes"]
        time_df[:,"nu"] .= loadinfo["nu"]
        time_df[:,"tot_samples"] .= result.sampled * 2
        loadinfo["out"] = "time"
        println("")
        println("-----------------------------------------------")
        println("")
        CSV.write(string("results/simulation/local/",simtypes[simnum],"-results/",savename(loadinfo,"csv",digits=1)),time_df)
    end
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
    
    psrf_df[:,"mean_xi"] .= mean(psrf.ξ[1:V])
    psrf_df[:,"max_xi"] .= max(psrf.ξ[1:V]...)
    psrf_df[:,"mean_gamma"] .= mean(psrf.γ[1:q])
    psrf_df[:,"max_gamma"] .= max(psrf.γ[1:q]...)

    @show psrf_df[1,"max_xi"];@show psrf_df[1,"max_gamma"]

    type = simtypes[simnum]
    
    saveinfo["out"] = "nodes"
    CSV.write(string("results/simulation/local/",type,"-results/",savename(saveinfo,"csv",digits=1)),output)

    saveinfo["out"] = "edges"
    CSV.write(string("results/simulation/local/",type,"-results/",savename(saveinfo,"csv",digits=1)),gam)

    saveinfo["out"] = "MSE"
    CSV.write(string("results/simulation/local/",type,"-results/",savename(saveinfo,"csv",digits=1)),mse_df)

    saveinfo["out"] = "mu"
    CSV.write(string("results/simulation/local/",type,"-results/",savename(saveinfo,"csv",digits=1)),mu_df)
    
    saveinfo["out"] = "psrf"
    CSV.write(string("results/simulation/local/",type,"-results/",savename(saveinfo,"csv",digits=1)),psrf_df)
    
end


main()
