using Test

include("../src/BayesNet.jl")


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
Z = symmetrize_matrices(X)
#show(stdout,"text/plain",Z)

η  = 1.01
ζ  = 1
ι  = 1
R  = 3
aΔ = 0
bΔ = 0
V = size(Z[1],1)
q = floor(Int,V*(V-1)/2)
n = size(Z,1)

@testset "InitTests" begin
    X, θ, D, πᵥ, Λ, Δ, ξ, M, u, μ, τ², γ, V_new = init_vars(Z, η, ζ, ι, R, aΔ, bΔ)
    #show(stdout, "text/plain",γ)
    #println("\n-\n")
    #show(stdout, "text/plain",γ[1])
    #println("\n-\n")
    #show(stdout, "text/plain", γ)
    @test (size(X,1),size(X[1],1)) == (n,q)
    @test size(D[1]) == (q,q)
    @test size(πᵥ[1]) == (R,3)
    @test size(Λ[1]) == (R,R)
    @test issubset(ξ[1],[0,1])
    @test size(M[1]) == (R,R)
    @test size(u[1]) == (R,V)
    @test size(γ[1]) == (q,)
    @test V == V_new
end
