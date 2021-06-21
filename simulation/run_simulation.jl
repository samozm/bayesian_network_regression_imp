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
    γ,MSE,ξ = sim_one_case(simnum,casenum,nburn,nsamp,jcon)

    if jcon
        ξ_in = DataFrame(CSV.File("juliacon/data/simulation$(simnum)_case$(casenum)_xis.csv"))
    else
        ξ_in = DataFrame(CSV.File("data/simulation/simulation$(simnum)_case$(casenum)_xis.csv"))

        step = Int(floor((nburn+nsamp)/100))
        ξ_prog = map(i->mean(ξ[1:i]),1:step:nburn+nsamp)
        #println("tst")
        #show(stdout,"text/plain",tst)
        #println("")
        savefig(plot(transpose(hcat(ξ_prog...)),legend=false),"plots/simulation/xi_converge_simulation$(simnum)_case$(casenum)")

        γ_prog = zeros(190,length(1:step:nburn+nsamp))
        for j in 1:190
            γ_prog[j,:] = map(i->mean(γ[j,1:i]),1:step:nburn+nsamp)
        end
        savefig(plot(transpose(γ_prog),legend=false),"plots/simulation/gamma_coverage_simulation$(simnum)_case$(casenum)")

    end

    output_results(γ[:,nburn+1:nburn+nsamp],MSE,mean(ξ[nburn+1:nburn+nsamp]),ξ_in,simnum,casenum,jcon)
end

function sim_one_case(simnum,casenum,nburn,nsamp,jcon,η=1.01,ζ=1,ι=1,R=5,aΔ=1,bΔ=1,ν=10)
    if jcon
        data_in = DataFrame(CSV.File("juliacon/data/simulation$(simnum)_case$(casenum).csv"))

        b_in = DataFrame(CSV.File("juliacon/data/simulation$(simnum)_case$(casenum)_bs.csv"))
        B₀ = convert(Array{Float64,1},b_in[!,:B])
    else
        data_in = DataFrame(CSV.File("data/simulation/simulation$(simnum)_case$(casenum).csv"))

        b_in = DataFrame(CSV.File("data/simulation/simulation$(simnum)_case$(casenum)_bs.csv"))
        B₀ = convert(Array{Float64,1},b_in[!,:B])
    end
    X = convert(Matrix,data_in[!,names(data_in,Not("y"))])
    y = data_in[:,:y]

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

function output_results(γ,MSE,ξ,ξ⁰,simnum,casenum,jcon)
    q = size(γ,1)
    V = convert(Int,(1 + sqrt(1 + 8*q))/2)

    plot_γ_sim(γ, "Gamma", "simulation$(simnum)_case$(casenum)",jcon)
    #ξ⁰[!,"Xi Probability"] = ξ
    #show(stdout,"text/plain",ξ⁰)
    #println("")
    show(stdout,"text/plain",DataFrame(MSE=MSE))
    println("")
    #println("Gamma")
    gam = DataFrame(create_upper_tri(vec(median.(γ)),V),:auto)
    #show(stdout,"text/plain",gam[!,names(gam,Not("x1"))])
    #println("")
    output = DataFrame(ξ⁰)
    #output[!,"MSE"] = MSE
    if jcon
        CSV.write("juliacon/results/simulation$(simnum)_case$(casenum).csv",output)
        CSV.write("juliacon/results/simulation$(simnum)_case$(casenum)_gammas.csv",gam)
    else
        CSV.write("results/simulation/simulation$(simnum)_case$(casenum).csv",output)
        CSV.write("results/simulation/simulation$(simnum)_case$(casenum)_gammas.csv",gam)
    end
end

main()
