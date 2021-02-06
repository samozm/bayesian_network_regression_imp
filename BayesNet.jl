using Random, Distributions, DataFrames, LinearAlgebra, StatsBase

"""
    sample_u(ξ, R, M)

Sample rows of the u matrix, either from MVN with mean 0 and covariance matrix M or a row of 0s

# Arguments
- `ξ` : ξ value sampled from Binomial distribution. Set to 0 to return a row of 0s, 1 to sample from MVN
- `R` : dimension of u vectors, length of return vector
- `M` : R×R covariance matrix for MVN samples 
"""
function sample_u(ξ, R, M)
    if (ξ == 1)
        return rand(MultivariateNormal(zeros(R), M))
    else
        return zeros(R)
    end
end

"""
    init_vars(X, η, ζ, ι, R, aΔ, bΔ)

Initialize all variables using prior distributions. Note, any value passed in a parameter marked as 'output parameter' will be ignored and overwritten.

# Arguments
- `X` : unweighted symmetric adjacency matrix to be used as a predictor
- `η` : hyperparameter used to sample from the Dirichlet distribution (r^η)
- `ζ` : hyperparameter used as the shape parameter in the gamma distribution used to sample θ
- `ι` : hyperparameter used as the scale parameter in the gamma distribution used to sample θ
- `R` : the dimensionality of the latent variables u, a hyperparameter
- `aΔ`: hyperparameter used as the a parameter in the beta distribution used to sample Δ. 
- `bΔ`: hyperparameter used as the b parameter in the beta distribution used to sample Δ. aΔ and bΔ values causing the Beta distribution to have mass concentrated closer to 0 will cause more zeros in ξ

# Returns
- `θ` : set to 1 draw from Gamma(ζ,ι)
- `s` : set to a V × V (where V × V is the dimensionality of X) matrix of draws from the Exponential(θ/2) distribution
- `πᵥ`: set to a R × 3 matrix of draws from the Dirichlet distribution, where the second and third columns are draws from Dirichlet(1) and the first are draws from Dirichlet(r^η)
- `Λ` : R × R diagonal matrix of λ values, which are sampled from [0,1,-1] with probabilities assigned from the rth row of πᵥ
- `Δ` : set to 1 draw from Beta(aΔ, bΔ) (or 1 if aΔ or bΔ are 0).
- `ξ` : set to V draws from Bernoulli(Δ)
- `M` : set to R × R matrix drawn from InverseWishart(V, I_R) (where I_R is the identity matrix with dimension R × R)
- `u` : the latent variables u, set to a V × R matrix where each row is sampled from MultivariateNormal(0,M) if ξ[r] is 1 or set to a row of 1s otherwise.
- `μ` : set to 1 (non-informative prior)
- `τ²`: set to the square of 1 sample from Uniform(0,1) (non-informative prior)
- `γ` : set to 1 draw from MultivariateNormal(uᵀΛu_upper, τ²*s), where uᵀΛu_upper is a vector of the terms in the upper triangle of uᵀΛu (does not include the diagonal)

See parameters marked 'output parameter'
"""
function init_vars(X, η, ζ, ι, R, aΔ, bΔ)
    # η must be greater than 1, if it's not set it to its default value of 1.01
    if (η <= 1)
        η = 1.01
        println("η value invalid, reset to default of 1.01")
    end
    V = size(X,1)
    q = V*(V-1)/2
    θ = rand(Gamma(ζ, ι))

    s = diagm(map(k -> rand(Exponential(θ/2)), 1:q))
    πᵥ = transpose(hcat(map(r -> rand(Dirichlet([r^η,1,1])),1:R)...))
    λ = map(r -> sample([0,1,-1], weights(πᵥ[r,:]),1)[1], 1:R)
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
    γ = rand(MultivariateNormal(uᵀΛu_upper, τ²*s))
    return (θ, s, πᵥ, Λ, Δ, ξ, M, u, μ, τ², γ)
end


function main()
    #data = DataFrame!(CSV.File())
    #X  = convert(Matrix, select(data, Not()))
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

    θ, s, πᵥ, Λ, Δ, ξ, M, u, μ, τ², γ = init_vars(X, η, ζ, ι, R, aΔ, bΔ)
end

main()

