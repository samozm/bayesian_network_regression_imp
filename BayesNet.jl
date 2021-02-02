using Random, Distributions, DataFrames, LinearAlgebra, StatsBase

function sample_u(ξ, R, M)
    if (ξ == 1)
        return rand(MultivariateNormal(zeros(R), M))
    else
        return Inf64 * ones(R)
    end
end

"""
    init_vars()

Initialize all variables using prior distributions

# Arguments

# Returns
"""
function init_vars(data, X, η, ζ, ι, R, aΔ, bΔ)
    if (η <= 1)
        η = 1.01
        println("η value invalid, reset to default of 1.01")
    end
    X = convert(Matrix, data)
    V = size(X,1)
    θ = rand(Gamma(ζ, ι))
    #TODO: there's probably a more efficient way to generate s
    s = convert(Matrix, reshape(map(k -> rand(Exponential(θ/2)), 1:V*V), V, V)) 
    π = map(r -> rand(Dirichlet([r^η,1,1])),1:R)
    λ = map(r -> sample([0,1,-1], weights(π[r]),1)[1], 1:R)
    Λ = diagm(λ)
    if (aΔ == 0 || bΔ == 0)
        Δ = 1
    else 
        Δ = rand(Beta(aΔ, bΔ))
    end
    ξ = map(keys -> rand(Binomial(1,Δ)), 1:V)
    M = rand(InverseWishart(V,Matrix(I,R,R)))
    u = hcat(map(k -> sample_u(ξ[k], R, M), 1:V)...)
    μ = 1
    τ²= rand(Uniform(0,1))^2
    uᵀΛu = transpose(u) * Λ * u
    uᵀΛu_upper = [uᵀΛu[i,j] for i = 1:size(uᵀΛu,1), j = 2:size(uᵀΛu,2) if i < j]
    γ = rand(MultivariateNormal(t))
end

function main()

end

main()

