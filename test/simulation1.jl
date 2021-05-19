using ArgParse,RCall,TickTock
include("../src/BayesNet.jl")
include("plot_output.jl")


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
    end
    return parse_args(args)
end

function main()
    parsed_CL_args = parse_CL_args()
    nburn = parsed_CL_args["nburn"]
    nsamp = parsed_CL_args["nsamp"]
    γ₁,MSE₁,ξ₁ = sim_one_case(nburn,nsamp)
    @rput nburn
    @rput nsamp
    R"source('test/run_guha_sim1.R');retlist <- guha_sim1(nburn,nsamp)"
    R"gamma <- retlist$gamma; MSE <- retlist$MSE; xis <- retlist$xis"
    γ₂ = @rget gamma
    MSE₂ = @rget MSE
    ξ₂ = @rget xis
    output_results(γ₁[:,nburn+1:nburn+nsamp],MSE₁,mean(ξ₁[nburn+1:nburn+nsamp]),γ₂,MSE₂,ξ₂)
end


function sim_one_case(nburn,nsamp)
    η  = 1.01
    ζ  = 1
    ι  = 1
    R  = 5
    aΔ = 1
    bΔ = 1
    ν = 12
    V = 20
    q = floor(Int,V*(V-1)/2)

    R"load('data/test/GuhaData.Rdata')"
    R"X <- simdata$Xmat; y <- simdata$y; b_in <- simdata$true.b"
    Z=@rget X
    y=@rget y
    B₀=@rget b_in
    println("Sim 1 Case 1")

    #region full run test
    tick()
    τ², u, ξ, γ, D, θ, Δ, M, μ, Λ, πᵥ = BayesNet(Z, y, R, nburn=nburn,nsamples=nsamp, V_in = 20, x_transform = false)
    tock()

    γ_n = hcat(γ...)


    γ_n2 = mean(γ[nburn+1:10:nburn+nsamp])
    γ₀ = B₀
    MSE = 0
    for i in 1:190
        MSE = MSE + (γ_n2[i] - γ₀[i])^2
    end
    MSE = MSE * (2/(V*(V-1)))

    return γ_n,MSE,ξ

end

function output_results(γ₁,MSE₁,ξ₁,γ₂,MSE₂,ξ₂)
    plot_γs(γ₁, γ₂, "Mine", "Guha", "sim1")
    show(stdout,"text/plain",DataFrame(My_Xi=ξ₁,Guha_Xi=ξ₂))
    println("")
    show(stdout,"text/plain",DataFrame(My_MSE=MSE₁,Guha_MSE=MSE₂))
    println("")
end

main()

#endregion


#region individual function tests
#=
R"outlist=readRDS('data/test/guha_out.Rdata')"
R"g_gamma=outlist$betamat; g_q=outlist$q; g_lambda=outlist$kappa; g_epsilon_k=outlist$epsilon_k; g_theta=outlist$theta; g_Lambda=outlist$Lambda; g_tau=outlist$tau; g_u=outlist$u; g_pi=outlist$pi; g_mu=outlist$mu; g_M=outlist$M; g_delta=outlist$delta; g_S=outlist$S; g_lambda=outlist$lambda"
g_mu=@rget g_mu
g_gamma=@rget g_gamma
g_tau=@rget g_tau
g_Lambda=@rget g_Lambda
g_u=@rget g_u
g_delta=@rget g_delta
g_epsilon=@rget g_epsilon_k
g_theta=@rget g_theta
g_lambda=@rget g_lambda
g_pi=@rget g_pi
g_M=@rget g_M
g_S=@rget g_S


println("mu")
mus = zeros(500)

for i in 1:500
    gam = upper_triangle(g_gamma[i + 32000,:,:])
    mus[i] = update_μ(Z, y, gam, g_tau[i + 32000,:][1], size(Z,1))
end

#show(stdout, "text/plain", DataFrame(n=collect(1:500), GuhaMu=g_mu[32001:32500], MeMu=mus))

GuhaSort=sort(g_mu[32001:32500])
MeSort=sort(mus)
lo = convert(Int64, round(500 * 0.10))
hi = convert(Int64, round(500 * 0.90))

show(stdout, "text/plain", DataFrame(GuhaLo=GuhaSort[lo], GuhaHi=GuhaSort[hi]))
println("")
show(stdout, "text/plain", DataFrame(MeLo=MeSort[lo], MeHi=MeSort[hi]))
println("")

println("tau")
taus = zeros(500)

for i in 1:500
    gam = upper_triangle(g_gamma[i + 32000,:,:])
    taus[i] = update_τ²(Z,y,g_mu[i + 32000],gam,diagm(g_lambda[i + 32000,:]),transpose(g_u[i + 32000]),diagm(upper_triangle(g_S[i + 32000,:,:])),20)
end

GuhaSort=sort(g_tau[32001:32500])
MeSort=sort(taus)

show(stdout, "text/plain", DataFrame(GuhaLo=GuhaSort[lo], GuhaHi=GuhaSort[hi]))
println("")
show(stdout, "text/plain", DataFrame(MeLo=MeSort[lo], MeHi=MeSort[hi]))
println("")

println("u,epsilon")
us = zeros(500)
epsilons = zeros(500)

println("gamma")
gammas = zeros(500,190)

for i in 1:500
    gammas[i,:] = update_γ(Z, diagm(upper_triangle(g_S[i + 32000,:,:])),diagm(g_lambda[i + 32000,:]),transpose(g_u[i + 32000]),g_mu[i + 32000],g_tau[i + 32000],70)
end

println(size(gammas))

g_gam_mn = zeros(190)
gam_mn = zeros(190)

for i in 1:100
    global g_gam_mn = g_gam_mn + upper_triangle(g_gamma[32000 + i,:,:])
    global gam_mn = gam_mn + gammas[i,:]
end
#show(stdout, "text/plain", DataFrame(guha=g_gam_mn .* (1/100),me=gam_mn .* (1/100)))
println("")

low = zeros(190)
high = zeros(190)
g_low = zeros(190)
g_high = zeros(190)
lw = convert(Int64, round(500 * 0.1))
hi = convert(Int64, round(500 * 0.9))

#γ_n = hcat(gammas...)
betasample = zeros(500,q)

for i in 1:500
    betasample[i,:] = upper_triangle(g_gamma[32000 + i,:,:])
end

println(size(betasample))

for i in 1:190
    srtd = sort(gammas[:,i])
    g_srtd = sort(betasample[:,i])
    low[i] = srtd[lw]
    high[i] = srtd[hi]
    g_low[i] = g_srtd[lw]
    g_high[i] = g_srtd[hi]
end

γ_df = DataFrame(n = collect(1:190),m_l = low, m_h = high, g_l = g_low, g_h = g_high)
show(stdout, "text/plain", γ_df)


println("D")
Ds = zeros(500,190,190)

for i in 1:500
    gam = upper_triangle(g_gamma[i + 32000,:,:])
    Ds[i,:,:] = update_D(gam,transpose(g_u[i + 32000]),diagm(g_lambda[i + 32000,:]),g_theta[i + 32000],g_tau[i + 32000],20)
end

g_gam_mn = zeros(190)
gam_mn = zeros(190)

for i in 1:100
    global g_gam_mn = g_gam_mn + upper_triangle(g_S[i + 32000,:,:])
    global gam_mn = gam_mn + diag(Ds[i,:,:])
end
show(stdout, "text/plain", DataFrame(guha=g_gam_mn .* (1/100),me=gam_mn .* (1/100)))
println("")

println("lambda")
lambdas = zeros(500,5,5)

for i in 1:500
    g_pi2 = g_pi[i + 32000,:,:]
    #g_pi2[:,1] = g_pi2[:,2]
    #g_pi2[:,2] = g_pi[i + 32000,:,1]
    gam = upper_triangle(g_gamma[i + 32000,:,:])
    lambdas[i,:,:] = update_Λ(g_pi2,R,diagm(g_lambda[i + 32000,:]),transpose(g_u[i + 32000]),diagm(upper_triangle(g_S[i + 32000,:,:])),g_tau[i + 32000],gam)
end

g_gam_mn = ones(1) .* 5
gam_mn = ones(1) .* 5

for i in 1:500
    #global g_gam_mn = g_gam_mn + g_lambda[i + 32000,:]
    append!(g_gam_mn, g_lambda[i+3200,:])
    append!(gam_mn,diag(lambdas[i,:,:]))
end
println(size(g_gam_mn))
println(size(gam_mn))

show(stdout, "text/plain", DataFrame(guha_0=count(x -> x==0, g_gam_mn) * (1/2500),
     guha_1=count(x -> x==1, g_gam_mn) * (1/2500),guha_neg1=count(x -> x==-1, g_gam_mn) * (1/2500),
     me_0=count(x -> x==0, gam_mn) * (1/2500),
     me_1=count(x -> x==1, gam_mn) * (1/2500), me_neg1=count(x -> x==-1, gam_mn) * (1/2500)))
println("")

println("theta")
thetas = zeros(500)

for i in 1:500
    thetas[i] = update_θ(ζ,ι,V,diagm(upper_triangle(g_S[i + 32000,:,:])))
end

show(stdout,"text/plain",size(g_theta[32001:32500]))
println("")

GuhaSort=sort(g_theta[32001:32500])
MeSort=sort(thetas)

show(stdout, "text/plain", DataFrame(GuhaLo=GuhaSort[lo], GuhaHi=GuhaSort[hi]))
println("")
show(stdout, "text/plain", DataFrame(MeLo=MeSort[lo], MeHi=MeSort[hi]))
println("")

println("delta")
deltas = zeros(500)

for i in 1:500
    deltas[i] = update_Δ(aΔ,bΔ,g_epsilon[i + 3200,:],20)
end

GuhaSort=sort(g_delta[32001:32500])
MeSort=sort(deltas)

show(stdout, "text/plain", DataFrame(GuhaLo=GuhaSort[lo], GuhaHi=GuhaSort[hi]))
println("")
show(stdout, "text/plain", DataFrame(MeLo=MeSort[lo], MeHi=MeSort[hi]))
println("")

println("M")
Ms = zeros(500,5,5)
uΛu = zeros(500,5,5)

for i in 1:500
    Ms[i,:,:] = update_M(transpose(g_u[32000 + i]),a_wish,20,g_epsilon[32010+i,:])
    #=
    if i == 10
        println("")
        show(stdout, "text/plain", Ms[i,:,:])
        println("")
    end
    =#
end
#show(stdout, "text/plain", mean(uΛu,dims=1))
println("")

g_m_mn = 0
m_mn = 0

for i in 1:500
    global g_m_mn = g_m_mn + sum(upper_triangle(g_M[i + 32000][:,:]))
    global m_mn = m_mn + sum(upper_triangle(Ms[i,:,:]))
end

show(stdout, "text/plain", DataFrame(guha=g_m_mn .* (1/100),me=m_mn .* (1/100)))
println("")

show(stdout,"text/plain",Ms[200,:,:])
println("")
show(stdout,"text/plain",g_M[32200,:,:])
println("")
println("pi")
pis = zeros(500,5,3)

for i in 1:500
    pis[i,:,:] = update_π(diagm(g_lambda[i + 32000,:]),η)
end

g_pi_mn = zeros(3)
pi_mn = zeros(3)

for i in 1:100
    for j in 1:3
        global g_pi_mn[j] = g_pi_mn[j] .+ sum(g_pi[i + 32000,:,j])
        global pi_mn[j] = pi_mn[j] .+ sum(pis[i,:,j])
    end
end
println("")
show(stdout, "text/plain", DataFrame(guha=g_pi_mn .* (1/100),me=pi_mn .* (1/100)))
println("")

println("u_ξ")
us = zeros(500,5,20)
epsilons = zeros(500,20)
mn = zeros(500,5,20)
cov=zeros(500,5,20)

for i in 1:500
    gam = upper_triangle(g_gamma[i + 32000,:,:])
    res = update_u_ξ(transpose(g_u[i + 32000]),gam,diagm(upper_triangle(g_S[i + 32000,:,:])),g_tau[i + 32000],g_delta[i + 32000],g_M[i + 32000],diagm(g_lambda[i + 32000,:]),20)
    us[i,:,:] = res[1]
    epsilons[i,:,:] = res[2]
    if i == 10 || i == 200 || i==400
        println(i)
        println("mu")
        show(stdout,"text/plain",res[3])
        println("cov")
        show(stdout,"text/plain",res[4])
    end
end


my_u = zeros(20,5)
gam_u = zeros(20,5)
my_ep = zeros(20)
gam_ep = zeros(20)
println(size(g_u[32000]))
for i in 1:20
    for j in 1:5
        my_u[i,j] = mean(us[:,j,i])
        #gam_u[i,j] = mean(g_u[32000:32500][j,i])
    end
    my_ep[i] = mean(epsilons[:,i])
    gam_ep[i] = mean(g_epsilon[32000:32500,i])
end
gam_u = mean(g_u[32000:32500])

#=
for i in [1,200,400]
    println(i)
    println("mean")
    show(stdout, "text/plain", mn[i,:,:])
    println("cov")
    show(stdout, "text/plain", cov[i,:,:])
end

println("my u")
show(stdout, "text/plain", my_u)
println("")
println("gamma u")
show(stdout, "text/plain", gam_u)
println("")
show(stdout, "text/plain", DataFrame(Me_ep=my_ep, Gamma_ep=gam_ep))
println("")
=#

#=

#println("y")
#show(stdout, "text/plain",y)
#println("")

#println("Xγ")
#show(stdout, "text/plain",Z*reshape(upper_triangle(g_gamma[32001,:,:]),190,1))
println("")
#println("X")
#show(stdout, "text/plain",Z[1,:,:])
println("")
#println("γ")
#show(stdout, "text/plain",upper_triangle(g_gamma[32001,:,:]))
println("")
=#

=#
#endregion
