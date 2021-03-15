using Test,TickTock,Suppressor

include("../src/BayesNet.jl")

#g_data = read_rda("data/GuhaData.Rdata")
R"g=load('data/GuhaData.Rdata')"
R"X <- simdata$Xmat; y <- simdata$y"
Z=@rget X
y=@rget y

function symmetrize_matrices(X)
    n = size(X,1)
    V = size(X[1],1)
    X_new = Array{Array{Int8,2},1}(undef,0)
    for i in 1:size(X,1)
        #show(stdout,"text/plain",X[i])
        B = convert(Matrix, reshape(X[i], 4, 4))
        #show(stdout,"text/plain",B)
        push!(X_new,Symmetric(B))
    end
    X = X_new
end

X = [[0, 1, 0, 1, 
     1, 0, 1, 1, 
     0, 1, 0, 0, 
     0, 1, 0, 0],
    [0, 1, 1, 1, 
     1, 0, 1, 1, 
     1, 1, 0, 0, 
     1, 1, 0, 0],
    [0, 0, 1, 0, 
     0, 0, 1, 1, 
     1, 1, 0, 0, 
     0, 1, 0, 0],
    [0, 0, 1, 0, 
     0, 0, 1, 1, 
     1, 1, 0, 1, 
     0, 1, 1, 0]]
#show(stdout,"text/plain",X)
#Z = symmetrize_matrices(X)
#show(stdout,"text/plain",Z)

#y = ones(size(Z[1],1))*12 + rand(Normal(0,2),size(Z[1],1))


η  = 1.01
ζ  = 1
ι  = 1
R  = 5
aΔ = 0
bΔ = 0
#V = size(Z,1)
V = 20
q = floor(Int,V*(V-1)/2)
n = size(Z,1)

#println("Z")
#show(stdout,"text/plain",Z)
#println(" ")


@testset "InitTests" begin
    X, θ, D, πᵥ, Λ, Δ, ξ, M, u, μ, τ², γ, V_new = init_vars(Z, η, ζ, ι, R, aΔ, bΔ,V,false)
    #show(stdout, "text/plain",γ)
    #println("\n-\n")
    #show(stdout, "text/plain",γ[1])
    #println("\n-\n")
    #show(stdout, "text/plain", γ)
    #show(stdout, "text/plain", Z)
    #show(stdout, "text/plain", X)
    @test size(X) == (n,q)
    @test size(D[1]) == (q,q)
    @test size(πᵥ[1]) == (R,3)
    @test size(Λ[1]) == (R,R)
    @test issubset(ξ[1],[0,1])
    @test size(M[1]) == (R,R)
    @test size(u[1]) == (R,V)
    @test size(γ[1]) == (q,)
    @test V == V_new
end
#=
nburn = 30000
nsamp = 20000
tick()
τ², u, ξ, γ, D, θ, Δ, M, μ, Λ, πᵥ = BayesNet(Z, y, R, nburn=nburn,nsamples=nsamp, V_in = 20, x_transform = false)
tock()

low = zeros(20)
high = zeros(20)
lw = round(20000 * 0.25)
hi = round(20000 * 0.75)
lw = convert(Int64, round(nsamp * 0.1))
hi = convert(Int64, round(nsamp * 0.9))

γ_n = hcat(γ...)

for i in 1:20
    srtd = sort(γ_n[i,nburn+1:nburn+nsamp])
    low[i] = srtd[lw]
    high[i] = srtd[hi]
end
=#
#show(stdout,"text/plain",DataFrame(n = collect(1:20),l = low, h = high))
#println("")




R"outlist=readRDS('data/guha_out.Rdata')"
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

println("g_M")
println(size(g_S))
println(size(g_S[1,:,:]))
#=
for i in 1:500
    print(size(upper_triangle(g_S[i + 32000,:,:])))
    print(size(g_delta))
    w=diagm(upper_triangle(g_S[i + 32000,:,:]))
    qqq=diagm(g_lambda[i + 32000,:])
    l=transpose(g_u[i + 32000])
    #print(typeof(g_M[i,:,:]))
    gam = upper_triangle(g_gamma[i + 32000,:,:])
    us[i],epsilons[i] = update_u_ξ(transpose(g_u[i + 32000]),gam,diagm(upper_triangle(g_S[i + 32000,:,:])),g_tau[i + 32000],g_delta[i + 32000],g_M[i + 32000],diagm(g_lambda[i + 32000,:]),20)
end

GuhaSort=sort(g_u[32001:32500])
MeSort=sort(us)

show(stdout, "text/plain", DataFrame(GuhaLo=GuhaSort[lo], GuhaHi=GuhaSort[hi]))
println("")
show(stdout, "text/plain", DataFrame(MeLo=MeSort[lo], MeHi=MeSort[hi]))
println("")
=#
GuhaSort=sort(g_epsilon[32001:32500])
MeSort=sort(epsilons)

show(stdout, "text/plain", DataFrame(GuhaLo=GuhaSort[lo], GuhaHi=GuhaSort[hi]))
println("")
show(stdout, "text/plain", DataFrame(MeLo=MeSort[lo], MeHi=MeSort[hi]))
println("")

println("gamma")
gammas = zeros(500,190)

for i in 1:500
    gammas[i,:,:] = update_γ(Z, diagm(upper_triangle(g_S[i + 32000,:,:])),diagm(g_lambda[i + 32000,:]),transpose(g_u[i + 32000]),g_mu[i + 32000],g_tau[i + 32000],70)
end

g_gam_mn = zeros(190)
gam_mn = zeros(190)

for i in 1:100
    global g_gam_mn = g_gam_mn + upper_triangle(g_gamma[32000 + i,:,:])
    global gam_mn = gam_mn + gammas[i,:]
end
show(stdout, "text/plain", DataFrame(guha=g_gam_mn .* (1/100),me=gam_mn .* (1/100)))
println("")


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

for i in 1:500
    Ms[i,:,:] = update_M(transpose(g_u[32000 + i]),diagm(g_lambda[i + 32000,:]),20)
end

g_m_mn = 0
m_mn = 0

for i in 1:100
    global g_m_mn = g_m_mn + sum(upper_triangle(g_M[i + 32000,:,:]))
    global m_mn = m_mn + sum(upper_triangle(Ms[i,:,:]))
end

show(stdout, "text/plain", DataFrame(guha=g_m_mn .* (1/100),me=m_mn .* (1/100)))
println("")

show(stdout,"text/plain",Ms[200,:,:])
show(stdout,"text/plain",g_M[32200,:,:])

##############
# START HERE #
##############
println("pi")
pis = zeros(500,5,5)

for i in 1:500
    pis[i,:] = update_π(diagm(g_lambda[i + 32000,:]),η)
end

GuhaSort=sort(upper_triangle(g_pi[32001:32500]))
MeSort=sort(pis)

show(stdout, "text/plain", DataFrame(GuhaLo=GuhaSort[lo], GuhaHi=GuhaSort[hi]))
println("")
show(stdout, "text/plain", DataFrame(MeLo=MeSort[lo], MeHi=MeSort[hi]))
println("")

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