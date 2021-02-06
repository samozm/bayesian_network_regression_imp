using Test

include("BayesNet.jl")

A = [0, 1, 0, 1, 
     1, 0, 1, 1, 
     0, 1, 0, 0, 
     0, 1, 0, 0]
B = convert(Matrix, reshape(A, 4, 4))
X = Symmetric(B)
η  = 1.01
ζ  = 1
ι  = 1
R  = 3
aΔ = 0
bΔ = 0
V = size(X,1)
q = V*(V-1)/2

@testset "InitTests" begin
    θ, s, πᵥ, Λ, Δ, ξ, M, u, μ, τ², γ = init_vars(X, η, ζ, ι, R, aΔ, bΔ)
    @test size(s) == (q,q)
    @test size(πᵥ) == (R,3)
    @test size(Λ) == (R,R)
    @test issubset(ξ,[0,1])
    @test size(M) == (R,R)
    @test size(u) == (R,V)
    @test size(γ) == (q,)
end