using Random, Distributions, DataFrames, LinearAlgebra, StatsBase, GenInvGaussian, RCall

#region custom sampling
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
    sample_rgig(b,k,θ)

Sample from the GeneralizedInverseGaussian distribution using RCall with p=1/2, b=b, a=a

# Arguments
- `a` : shape and scale parameter a, sometimes also called ψ
- `b` : shape and scale parameter b, sometimes also called χ

# Returns
one sample from the GIG distribution with p=1/2, b=b, a=a
"""
function sample_rgig(a,b)
    @rput b
    @rput a
    R"r=rgig(n=1, lambda=1/2,chi=a,psi=b)"
    a=@rget r
end
#endregion

#region beginning of implementing sampling from GIG (unfinished)
function sample_GIG(λ,Ψ,χ)
    β = sqrt(Ψ*χ)
    if(β > 1)
        println("Uh oh, bad parameters for sample GIG")
        return
    elseif(β < min(1/2, (2/3)*sqrt(1 - λ))) 
        return sample_GIG_bigβ(λ,β)
    else
        return sample_GIG_smallβ(λ,β)
    end
end

function sample_GIG_bigβ(λ,β)
    m = β/((1 - λ) + sqrt((1 - λ)^2 + β^2))
    x⁺= ((1 + λ) + sqrt((1 + λ)^2 + β^2))/β
    v⁺= sqrt(g(m,λ,β))
    u⁺= x⁺ * sqrt(g(x⁺,λ,β))
    V = 0
    while true
        
    end
end

function sample_GIG_smallβ(λ,β)
    
end

function g(x,λ,β)
    x^(λ-1)*exp(-(β/2)*(x + (1/x)))
end

#endregion

#region helper functions
"""
    upper_triangle(matrix)

return the upper triangle (without the diagonal) of the matrix as a vector

# Arguments
- `matrix`: matrix of which to capture the upper triangle

# Returns
vector of upper triangluar section of `matrix`
"""
function upper_triangle(matrix)
    [matrix[i,j] for i = 1:size(matrix,1), j = 2:size(matrix,2) if i < j]
end
#endregion

"""
    init_vars(X, η, ζ, ι, R, aΔ, bΔ)

    Initialize all variables using prior distributions. Note, any value passed in a parameter marked as 'output parameter' will be ignored and overwritten.

    # Arguments
    - `X` : vector of unweighted symmetric adjacency matrices to be used as predictors. each element of the array should be 1 matrix
    - `η` : hyperparameter used to sample from the Dirichlet distribution (r^η)
    - `ζ` : hyperparameter used as the shape parameter in the gamma distribution used to sample θ
    - `ι` : hyperparameter used as the scale parameter in the gamma distribution used to sample θ
    - `R` : the dimensionality of the latent variables u, a hyperparameter
    - `aΔ`: hyperparameter used as the a parameter in the beta distribution used to sample Δ. 
    - `bΔ`: hyperparameter used as the b parameter in the beta distribution used to sample Δ. aΔ and bΔ values causing the Beta distribution to have mass concentrated closer to 0 will cause more zeros in ξ

    # Returns
    - `X` : matrix of re-ordered predictors. one row per sample, V*(V-1) columns 
    - `θ` : set to 1 draw from Gamma(ζ,ι)
    - `D` : set to a diagonal matrix of V(V-1) draws from the Exponential(θ/2) distribution
    - `πᵥ`: set to a R × 3 matrix of draws from the Dirichlet distribution, where the second and third columns are draws from Dirichlet(1) and the first are draws from Dirichlet(r^η)
    - `Λ` : R × R diagonal matrix of λ values, which are sampled from [0,1,-1] with probabilities assigned from the rth row of πᵥ
    - `Δ` : set to 1 draw from Beta(aΔ, bΔ) (or 1 if aΔ or bΔ are 0).
    - `ξ` : set to V draws from Bernoulli(Δ)
    - `M` : set to R × R matrix drawn from InverseWishart(V, I_R) (where I_R is the identity matrix with dimension R × R)
    - `u` : the latent variables u, set to a V × R matrix where each row is sampled from MultivariateNormal(0,M) if ξ[r] is 1 or set to a row of 1s otherwise.
    - `μ` : set to 1 (non-informative prior)
    - `τ²`: set to the square of 1 sample from Uniform(0,1) (non-informative prior)
    - `γ` : set to 1 draw from MultivariateNormal(uᵀΛu_upper, τ²*s), where uᵀΛu_upper is a vector of the terms in the upper triangle of uᵀΛu (does not include the diagonal)
    - `V` : dimension of original symmetric adjacency matrices
"""
function init_vars(X, η, ζ, ι, R, aΔ, bΔ)
    # η must be greater than 1, if it's not set it to its default value of 1.01
    if (η <= 1)
        η = 1.01
        println("η value invalid, reset to default of 1.01")
    end

    V = size(X[1],1)
    q = floor(Int,V*(V-1)/2)

    X_new = Matrix{Float64}(undef, size(X,1), q)
    X_new = map(k -> upper_triangle(X[k]), 1:size(X,1))
    
    θ = rand(Gamma(ζ, ι))

    D = diagm(map(k -> rand(Exponential(θ/2)), 1:q))
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
    uᵀΛu_upper = upper_triangle(uᵀΛu)
    γ = rand(MultivariateNormal(uᵀΛu_upper, τ²*D))
    return (X_new, θ, D, πᵥ, Λ, Δ, ξ, M, u, μ, τ², γ, V)
end

#region update variables
"""
    update_μ(y, X, γ, τ², n)

Sample the next μ value from the normal distribution with mean 1ᵀ(y - Xγ)/n and variance τ²/n

# Arguments
- `y` : response values 
- `X` : 2 dimensional array of predictor values, 1 row per sample (upper triangle of original X)
- `γ` : vector of regression parameters 
- `τ²`: overall variance parameter
- `n` : number of samples (length of y)

# Returns
new value of μ
"""
function update_μ(y, X, γ, τ², n)
    μₘ = (ones(1,n) * (y - X*γ)) / n
    σₘ² = τ²/n
    μ = rand(Normal(μₘ,σₘ²))
end

"""
    update_γ(X, D, W, μ, τ²)

Sample the next γ value from the normal distribution, decomposed as described in Guha & Rodriguez 2018

# Arguments
- `X` : 2 dimensional array of predictor values, 1 row per sample (upper triangle of original X)
- `D` : diagonal matrix of s values
- `W` : vector of uΛu values (upper triangle of uΛu), measures the effect of the relationship between nodes on the response
- `μ` : overall mean value for the relationship
- `τ²`: overall variance parameter 
- `n` : number of samples 

# Returns
new value of γ
"""
function update_γ(X, D, W, μ, τ², n)
    Δᵧ₁ = rand(MultivariateNormal(zeros(n,1), τ²*D))
    Δᵧ₂ = rand(MultivariateNormal(zeros(n,1), I(n)))
    Δᵧ₃ = (X/sqrt(τ²))*Δᵧ₁ + Δᵧ₂
    γw = Δᵧ₁ + (τ²*D)*(transpose(X)/sqrt(τ²))*inv(X*D*transpose(X) + I(n)) * (((y - μ*ones(n,1) - X*W)/sqrt(τ²)) - Δᵧ₃)
    γ = γw + W
end

"""
    update_τ²(X, y, μ, γ, W, D, V) 

Sample the next τ² value from the InverseGaussian distribution with mean n/2 + V(V-1)/4 and variance ((y - μ1 - Xγ)ᵀ(y - μ1 - Xγ) + (γ - W)ᵀD⁻¹(γ - W)

# Arguments
- `X` : 2 dimensional array of predictor values, 1 row per sample (upper triangle of original X)
- `y` : vector of response values
- `μ` : overall mean value for the relationship
- `γ` : vector of regression parameters 
- `W` : vector of uΛu values (upper triangle of uΛu), measures the effect of the relationship between nodes on the response
- `D` : diagonal matrix of s values
- `V` : dimension of original symmetric adjacency matrices

# Returns
new value of τ²
"""
function update_τ²(X, y, μ, γ, W, D, V)
    n  = size(y)
    μₜ  = (n/2) + (V*(V-1)/4)
    yμ1Xγ = (y - μ*ones(n,1) - X*γ)
    γW = (γ - W)
    σₜ² = (transpose(yμ1Xγ)*yμ1Xγ + transpose(γW)*inv(D)*γW)/2
    τ² = rand(InverseGaussian(μₜ, σₜ²))
end

"""
    update_D(γ, u, Λ, τ², V)

Sample the next D value from the GeneralizedInverseGaussian distribution with a = 1/2, b=((γ - uᵀΛu)^2)/τ², p=θ

# Arguments
- `γ` : vector of regression parameters 
- `u` : the latent variables u
- `Λ` : R × R diagonal matrix of λ values
- `τ²`: overall variance parameter 
- `V` : dimension of original symmetric adjacency matrices

# Returns
new value of D
"""
function update_D(γ, u, Λ, τ², V)
    q = floor(Int,V*(V-1)/2)
    uᵀΛu = transpose(u) * Λ * u
    uᵀΛu_upper = upper_triangle(uᵀΛu)
    a = (γ - uᵀΛu_upper).^2 / τ²
    D = diagm(map(k -> sample_rgig(a[k],θ), 1:q))
end

function update_θ()
    
end

function update_u()
    
end

function update_ξ()
    
end

function update_Δ()
    
end

function update_M()
    
end

function update_λ()
    
end

function update_π()

end
#endregion

function main()
    #data = DataFrame!(CSV.File())
    #X  = convert(Matrix, select(data, Not()))
    #A = [0, 1, 0, 1, 
    #    1, 0, 1, 1, 
    #    0, 1, 0, 0, 
    #    0, 1, 0, 0]
    #B = convert(Matrix, reshape(A, 4, 4))
    #X = Symmetric(B)
    η  = 1.01
    ζ  = 1
    ι  = 1
    R  = 3
    aΔ = 0
    bΔ = 0

    θ, s, πᵥ, Λ, Δ, ξ, M, u, μ, τ², γ = init_vars(X, η, ζ, ι, R, aΔ, bΔ)
end

#main()

