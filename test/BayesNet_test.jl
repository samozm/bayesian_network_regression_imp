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

nburn = 30000
nsamp = 20000
tick()
#τ², u, ξ, γ, D, θ, Δ, M, μ, Λ, πᵥ = BayesNet(Z, y, R, nburn=nburn,nsamples=nsamp, V_in = 20, x_transform = false)
tock()

low = zeros(20)
high = zeros(20)
#lw = round(20000 * 0.25)
#hi = round(20000 * 0.75)
lw = convert(Int64, round(nsamp * 0.1))
hi = convert(Int64, round(nsamp * 0.9))

#γ_n = hcat(γ...)

#for i in 1:20
#    srtd = sort(γ_n[i,nburn+1:nburn+nsamp])
#    low[i] = srtd[lw]
#    high[i] = srtd[hi]
#end

#show(stdout,"text/plain",DataFrame(n = collect(1:20),l = low, h = high))
println("")


R"outlist=readRDS('data/guha_out.Rdata')"
R"g_gamma=outlist$betamat; g_q=outlist$q; g_lambda=outlist$kappa; g_epsilon_k=outlist$epsilon_k; g_theta=outlist$theta; g_Lambda=outlist$Lambda; g_tau=outlist$tau; g_u=outlist$u; g_pi=outlist$pi; g_mu=outlist$mu; g_M=outlist$M; g_delta=outlist$delta; g_S=outlist$S"
g_mu=@rget g_mu
g_gamma=@rget g_gamma
g_tau=@rget g_tau

mus = zeros(500)

for i in 1:500
    gam = upper_triangle(g_gamma[i + 32000,:,:])
    mus[i] = update_μ(Z, y, gam, g_tau[i + 32000,:][1], size(Z,1))
end

show(stdout, "text/plain", DataFrame(n=collect(1:500), GuhaMu=g_mu[32001:32500], MeMu=mus))

GuhaSort=sort(g_mu[32001:32500])
MeSort=sort(mus)
lo = convert(Int64, round(500 * 0.10))
hi = convert(Int64, round(500 * 0.90))

show(stdout, "text/plain", DataFrame(GuhaLo=GuhaSort[lo], GuhaHi=GuhaSort[hi]))
show(stdout, "text/plain", DataFrame(GuhaMu=MeSort[lo], MeMu=MeSort[hi]))

