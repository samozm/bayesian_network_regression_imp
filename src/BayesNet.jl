using Random, DataFrames, LinearAlgebra, StatsBase, RCall, InvertedIndices
using Distributions
R"library(GIGrvg)"

#include("../../../Distributions.jl/src/Distributions.jl")

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
    sample_rgig(a,b)

Sample from the GeneralizedInverseGaussian distribution using RCall with p=1/2, b=b, a=a

# Arguments
- `a` : shape and scale parameter a, sometimes also called ψ
- `b` : shape and scale parameter b, sometimes also called χ

# Returns
one sample from the GIG distribution with p=1/2, b=b, a=a
"""
function sample_rgig(a,b)
    #TODO: confirm documentation/parameters are correct
    #@rput b
    #@rput a
    #R"r=rgig(n=1, lambda=1/2,chi=b,psi=a)"
    #l=@rget r
    
    #print("a: ")
    #println(a)
    #print("b: ")
    #println(b)
    l = rand(GeneralizedInverseGaussian(a,b,1/2))
    #println("end GIG")
    return l
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

"""
    create_upper_tri(vec,V)

create an upper triangluar matrix from a vector of the form [12, ... 1V,23,...(V-1)V] 
to the form [0 12 13  ... 1V]
            [0 0  23  ... 2V]
            [.............]
            [0 0  0...(V-1)V]
# Arguments
- `vec`: vector containing values to put into the upper triangluar matrix
- `V`  : dimension of output matrix

# Returns
upper triangluar matrix containing values of `vec` 
"""
function create_upper_tri(vec,V)
    mat = zeros(V,V)
    vec2 = deepcopy(vec)
    for k = 1:V
        for l = k+1:V
            mat[k,l] = popfirst!(vec2)
        end
    end
    return mat
end


function sample_π_dirichlet(r,η,λ)
    wts = [r^η,1,1]
    if(λ[r] == 0)
        wts[1] = r^η + 1
    elseif(λ[r] == 0)
        wts[2] = 2
    else
        wts[3] = 2
    end
    rand(Dirichlet(wts))
end
#endregion

"""
    init_vars(X, η, ζ, ι, R, aΔ, bΔ, V_in, x_transform)

    Initialize all variables using prior distributions. Note, any value passed in a parameter marked as 'output parameter' will be ignored and overwritten.

    # Arguments
    - `X` : vector of unweighted symmetric adjacency matrices to be used as predictors. each element of the array should be 1 matrix
    - `η` : hyperparameter used to sample from the Dirichlet distribution (r^η)
    - `ζ` : hyperparameter used as the shape parameter in the gamma distribution used to sample θ
    - `ι` : hyperparameter used as the scale parameter in the gamma distribution used to sample θ
    - `R` : the dimensionality of the latent variables u, a hyperparameter
    - `aΔ`: hyperparameter used as the a parameter in the beta distribution used to sample Δ. 
    - `bΔ`: hyperparameter used as the b parameter in the beta distribution used to sample Δ. aΔ and bΔ values causing the Beta distribution to have mass concentrated closer to 0 will cause more zeros in ξ
    - `V_in`: Value of V, the number of nodes in the original X matrix. Only used when x_transform is false
    - `x_transform`: boolean, set to false if X has been pre-transformed into one row per sample. True by default.

    # Returns
    - `X` : matrix of re-ordered predictors. one row per sample, V*(V-1)/2 columns 
    - `θ` : set to 1 draw from Gamma(ζ,ι)
    - `D` : set to a diagonal matrix of V(V-1)/2 draws from the Exponential(θ/2) distribution
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
function init_vars(X, η, ζ, ι, R, aΔ, bΔ, V_in=NaN, x_transform=true)
    # η must be greater than 1, if it's not set it to its default value of 1.01
    if (η <= 1)
        η = 1.01
        println("η value invalid, reset to default of 1.01")
    end

    if x_transform
        V = size(X[1],1)
    else
        V = V_in
    end 
    q = floor(Int,V*(V-1)/2)

    if x_transform
        X_new = Matrix{Float64}(undef, size(X,1), q)
        X_new = transpose(hcat(map(k -> upper_triangle(X[k]), 1:size(X,1))...))
    else
        X_new = X
    end

    θ = rand(Gamma(ζ, ι))

    S = map(k -> rand(Exponential(θ/2)), 1:q)
    D = diagm(S)
    πᵥ = transpose(hcat(map(r -> rand(Dirichlet([r^η,1,1])),1:R)...))
    λ = map(r -> sample([0,1,-1], weights(πᵥ[r,:]),1)[1], 1:R)
    Λ = diagm(λ)
    if (aΔ == 0 || bΔ == 0)
        Δ = 0.5
    else 
        Δ = rand(Beta(aΔ, bΔ))
    end
    #println("V")
    #println(V)
    #println("R")
    #println(R)

    ξ = map(keys -> rand(Binomial(1,Δ)), 1:V)
    M = rand(InverseWishart(V,Matrix(I,R,R)))
    u = hcat(map(k -> sample_u(ξ[k], R, M), 1:V)...)
    μ = 1.0
    τ²= rand(Uniform(0,1))^2
    uᵀΛu = transpose(u) * Λ * u
    uᵀΛu_upper = upper_triangle(uᵀΛu)
    γ = rand(MultivariateNormal(uᵀΛu_upper, τ²*D))
    return (X_new, [θ], [D], [πᵥ], [Λ], [Δ], [ξ], [M], [u], [μ], [τ²], [γ], V)
end

#region update variables
"""
    update_μ(y, X, γ, τ², n)

Sample the next μ value from the normal distribution with mean 1ᵀ(y - Xγ)/n and variance τ²/n

# Arguments
- `X` : 2 dimensional array of predictor values, 1 row per sample (upper triangle of original X)
- `y` : response values 
- `γ` : vector of regression parameters 
- `τ²`: overall variance parameter
- `n` : number of samples (length of y)

# Returns
new value of μ
"""
function update_μ(X, y, γ, τ², n)
    μₘ = (ones(1,n) * (y - X*γ)) / n
    σₘ² = τ²/n
    μ = rand(Normal(μₘ[1],σₘ²))
end

"""
    update_γ(X, D, Λ, u, μ, τ²)

Sample the next γ value from the normal distribution, decomposed as described in Guha & Rodriguez 2018

# Arguments
- `X` : 2 dimensional array of predictor values, 1 row per sample (upper triangle of original X)
- `D` : diagonal matrix of s values
- `Λ` : R × R diagonal matrix of λ values
- `u` : the latent variables u
- `μ` : overall mean value for the relationship
- `τ²`: overall variance parameter 
- `n` : number of samples 

# Returns
new value of γ
"""
function update_γ(X, D, Λ, u, μ, τ², n)
    uᵀΛu = transpose(u) * Λ * u
    W = upper_triangle(uᵀΛu)
    q = size(D,1)

    Δᵧ₁ = rand(MultivariateNormal(zeros(q), τ²*D))
    Δᵧ₂ = rand(MultivariateNormal(zeros(n), I(n)))
    Δᵧ₃ = (X/sqrt(τ²))*Δᵧ₁ + Δᵧ₂
    γw = Δᵧ₁ + (τ²*D)*(transpose(X)/sqrt(τ²))*inv(X*D*transpose(X) + I(n)) * (((y - μ*ones(n,1) - X*W)/sqrt(τ²)) - Δᵧ₃)
    γ = γw + W
end

"""
    update_τ²(X, y, μ, γ, Λ, u, D, V) 

Sample the next τ² value from the InverseGaussian distribution with mean n/2 + V(V-1)/4 and variance ((y - μ1 - Xγ)ᵀ(y - μ1 - Xγ) + (γ - W)ᵀD⁻¹(γ - W)

# Arguments
- `X` : 2 dimensional array of predictor values, 1 row per sample (upper triangle of original X)
- `y` : vector of response values
- `μ` : overall mean value for the relationship
- `γ` : vector of regression parameters 
- `Λ` : R × R diagonal matrix of λ values
- `u` : the latent variables u
- `D` : diagonal matrix of s values
- `V` : dimension of original symmetric adjacency matrices

# Returns
new value of τ²
"""
function update_τ²(X, y, μ, γ, Λ, u, D, V)
    uᵀΛu = transpose(u) * Λ * u
    W = upper_triangle(uᵀΛu)
    n  = size(y,1)
    
    #TODO: better variable names, not so much reassignment
    μₜ  = (n/2) + (V*(V-1)/4)

    #println("y")
    #show(stdout, "text/plain", y)
    #println("")
    #println("μ")
    #show(stdout, "text/plain", μ)
    #println("")
    #println("X*γ")
    #show(stdout, "text/plain", X*γ)
    #println("")

    #TODO: why is every yμ1Xγ pretty much the same?
    yμ1Xγ = (y .- μ*ones(n,1) .- X*γ)
    #println("yμ1Xγ")
    #show(stdout, "text/plain", yμ1Xγ)
    #println("")
    γW = (γ .- W)
    yμ1Xγᵀyμ1Xγ = transpose(yμ1Xγ)*yμ1Xγ
    γWᵀγW = transpose(γW)*inv(D)*γW

    #println("yμ1Xγᵀyμ1Xγ")
    #show(stdout, "text/plain", yμ1Xγᵀyμ1Xγ)
    #println("")

    σₜ² = (yμ1Xγᵀyμ1Xγ[1] + γWᵀγW[1])/2
    τ² = rand(InverseGamma(μₜ, σₜ²))
end

"""
    update_D(γ, u, Λ, θ, τ², V)

Sample the next D value from the GeneralizedInverseGaussian distribution with p = 1/2, a=((γ - uᵀΛu)^2)/τ², b=θ

# Arguments
- `γ` : vector of regression parameters 
- `u` : the latent variables u
- `Λ` : R × R diagonal matrix of λ values
- `θ` : b parameter for the GeneralizedInverseGaussian distribution
- `τ²`: overall variance parameter 
- `V` : dimension of original symmetric adjacency matrices

# Returns
new value of D
"""
function update_D(γ, u, Λ, θ, τ², V)
    q = floor(Int,V*(V-1)/2)
    uᵀΛu = transpose(u) * Λ * u
    uᵀΛu_upper = upper_triangle(uᵀΛu)
    a = (γ - uᵀΛu_upper).^2 / τ²
    # TODO: update this to use new sampler from distributions
    D = diagm(map(k -> sample_rgig(θ,a[k]), 1:q))
end

"""
    update_θ(ζ, ι, V, D)

Sample the next θ value from the Gamma distribution with a = ζ + V(V-1)/2 and b = ι + ∑(s[k,l]/2)

# Arguments
- `ζ` : hyperparameter, used to construct `a` parameter
- `ι` : hyperparameter, used to construct `b` parameter
- `V` : dimension of original symmetric adjacency matrices
- `D` : diagonal matrix of s values
- `τ²`: overall variance parameter 
# Returns
new value of θ
"""
function update_θ(ζ, ι, V, D)
    a = ζ + (V*(V-1))/2
    b = ι + sum(D)/2
    θ = rand(Gamma(a,b))
end

"""
    update_u_ξ(u, γ, D, τ², Δ, M, Λ, V)

Sample the next u and ξ values 

# Arguments
- `u` : the current latent variables u
- `γ` : vector of regression parameters 
- `D` : diagonal matrix of s values
- `τ²`: overall variance parameter
- `Δ` : set to 1 draw from Beta(aΔ, bΔ) (or 1 if aΔ or bΔ are 0).
- `M` : set to R × R matrix drawn from InverseWishart(V, I_R) (where I_R is the identity matrix with dimension R × R)
- `Λ` : R × R diagonal matrix of λ values
- `V` : dimension of original symmetric adjacency matrices

# Returns
A touple with the new values of u and ξ
"""
function update_u_ξ(u, γ, D, τ², Δ, M, Λ, V)
    q = V*(V-1)
    w_top = zeros(V)
    u_new = zeros(size(u)...)
    ξ = zeros(V)
    for k in 1:V
        #println("size(u)")
        #println(size(u))
        #println("size(Λ)")
        #println(size(Λ))
        U = transpose(u[:,Not(k)])*Λ
        s = create_upper_tri(diag(D),V)
        Γ = create_upper_tri(γ, V)
        if k == 1
            γk=vcat(Γ[1,2:V])
            H = diagm(vcat(s[1,2:V]))
        elseif k == V
            γk=vcat(Γ[1:V-1,V])
            H = diagm(vcat(s[1:V-1,V]))
        else
            H = diagm(vcat(s[1:k-1,k],s[k,k+1:V]))
            γk= vcat(Γ[1:k-1,k],Γ[k,k+1:V])
        end
        #println("D")
        #show(stdout, "text/plain", D[1:5,1:5])
        #println("")
        #println("s")
        #show(stdout, "text/plain", s[1:5,1:5])
        #println("")
        #println("H")
        #show(stdout, "text/plain", H[1:5,1:5])
        #println("")
        Σ = inv(((transpose(U)*inv(H)*U)/τ²) + inv(M))
        m = (Σ*transpose(U)*inv(H)*γk)/τ²
        #println("τ²")
        #show(stdout, "text/plain", τ²)
        #println("")
        #println("H")
        #show(stdout, "text/plain", H)
        #println("")
        mvn_a = MultivariateNormal(zeros(size(H,1)),τ²*H)
        mvn_b_Σ = round.(τ²*H + U*M*transpose(U), digits=10)
        mvn_b = MultivariateNormal(zeros(size(H,1)),mvn_b_Σ)
        w_top = (1-Δ) * pdf(mvn_a,γk)
        w_bot = w_top + Δ*pdf(mvn_b,γk)
        if (w_bot == 0)
            println("τ²")
            show(stdout, "text/plain", τ²)
            println("")
            println("H")
            show(stdout, "text/plain", H)
            println("")
            println("w")
            show(stdout, "text/plain", w)
            println("")
        end
        w = w_top / w_bot
        mvn_f = MultivariateNormal(m,round.(Σ,digits=10))
        
        
        #TODO: sometimes w is NaN
        ξ[k] = update_ξ(w)
        # the paper says the first term is (1-w) but their code uses 1-ξ. again i think this makes more sense
        # that this term would essentially be an indicator rather than a weight
        u_new[:,k] = (1-ξ[k]).*rand(mvn_f)
    end
    return u_new,ξ
end

"""
    update_ξ(w)

Sample the next ξ value from the Bernoulli distribution with parameter 1-w

# Arguments
- `w` : parameter to use for sampling, probability that 0 is drawn

# Returns
the new value of ξ
"""
function update_ξ(w)
    rand(Bernoulli(1 - w))
end

"""
    update_Δ(aΔ, bΔ, ξ, V)

Sample the next Δ value from the Beta distribution with parameters a = aΔ + ∑ξ and b = bΔ + V - ∑ξ

# Arguments
- `aΔ`: hyperparameter used as part of the a parameter in the beta distribution used to sample Δ. 
- `bΔ`: hyperparameter used as part of the b parameter in the beta distribution used to sample Δ.
- `ξ` : vector of ξ values, 0 or 1
- `V` : dimension of original symmetric adjacency matrices

# Returns
the new value of Δ
"""
function update_Δ(aΔ, bΔ, ξ, V)
    a = aΔ + sum(ξ)
    #TODO: ensure b is > 0
    b = bΔ + V - sum(ξ)
    if (a == 0 || b == 0)
        Δ = 0.5
    else 
        Δ = rand(Beta(a, b))
    end
    return Δ
end


"""
    update_M(u,Λ,V)

Sample the next M value from the InverseWishart distribution with df = V + # of nonzero columns in u and Ψ = I + ∑ uΛuᵀ

# Arguments
- `u` : R × V matrix of latent variables
- `Λ` : R × R diagonal matrix of λ values
- `V` : dimension of original symmetric adjacency matrices

# Returns
the new value of M
"""
function update_M(u,Λ,V)
    R = size(u,1)
    uΛu = 0
    num_nonzero = 0
    for i in 1:V
        # TODO: in the paper this is uᵀΛu but in the code it's just uuᵀ
        # should we email them and ask?
        # once the code works better I can play with the difference a bit more and see what works better
        uΛu = uΛu .+ u[:,i]*transpose(u[:,i])#*Λ*u[:,i]
        if u[:,i] != zeros(size(u,1))
            num_nonzero = num_nonzero + 1
        end
    end
    Ψ = I(R) .+ uΛu
    df = V + num_nonzero
    M = rand(InverseWishart(df,round.(Ψ, digits=5)))
    return M
end

"""
    update_Λ(πᵥ, R, λ, u, τ², D)

Sample the next values of λ from [0,1,-1] with probabilities determined from a normal mixture

# Arguments
- `πᵥ`: 3 column vectors of dimension R used to weight normal mixture for probability values
- `R` : the dimensionality of the latent variables u
- `Λ` : R × R diagonal matrix of λ values
- `u` : R × V matrix of latent variables
- `D` : diagonal matrix of s values
- `τ²`: overall variance parameter
- `γ` : vector of regression parameters 

# Returns
new value of Λ
"""
function update_Λ(πᵥ, R, Λ, u, D, τ², γ)
    λ = diag(Λ)
    λ_new = zeros(size(Λ,1))
    for r in 1:R
        Λ₋₁= Λ
        Λ₋₁[r,r] = -1
        Λ₀ = Λ
        Λ₀[r,r] = 0
        Λ₁ = Λ
        Λ₀[r,r] = 1
        # TODO: confirm it's correct to have the transpose u first
        # is this actually asking for the pdfs here?
        W₋₁= upper_triangle(transpose(u) * Λ₋₁ * u)
        W₀ = upper_triangle(transpose(u) * Λ₀ * u)
        W₁ = upper_triangle(transpose(u) * Λ₁ * u)
        n₀ = pdf(MultivariateNormal(W₀, τ² * D),γ)
        n₁ = pdf(MultivariateNormal(W₁, τ² * D),γ)
        n₋₁= pdf(MultivariateNormal(W₋₁, τ² * D),γ)
        p_bot = πᵥ[r,1] * n₀ + πᵥ[r,2] * n₁ + πᵥ[r,3] * n₋₁
        p1 = πᵥ[r,1] * n₀ / p_bot
        p2 = πᵥ[r,2] * n₁ / p_bot
        p3 = 1 - p1 - p2
        λ_new[r] = sample([0,1,-1],weights([p1,p2,p3]))
    end
    return diagm(λ_new)
end

"""
    update_π(λ,η)

Sample the new values of πᵥ from the Dirichlet distribution with parameters [#{r: λᵣ = 0} + r^η, 1 + #{r: λᵣ= 1}, 1 + #{r: λᵣ = -1 }]

# Arguments
- `Λ` : R × R diagonal matrix of λ values
- `η` : hyperparameter used for sampling the 0 value (r^η)

# Returns
new value of πᵥ
"""
function update_π(Λ,η)
    λ = diag(Λ)
    π_new = transpose(hcat(map(r -> sample_π_dirichlet(r,η,λ),1:R)...))
end
#endregion

"""
    GibbsSample(X, y, θ, D, πᵥ, Λ, Δ, ξ, M, u, μ, γ, V, η, ζ, ι, R, aΔ, bΔ)

Take one GibbsSample

# Arguments


# Returns
A tuple of the new values, (τ²_n, u_n, ξ_n, γ_n, D_n, θ_n, Δ_n, M_n, μ_n, Λ_n, πᵥ_n)
"""
function GibbsSample(X, y, θ, D, πᵥ, Λ, Δ, ξ, M, u, μ, γ, V, η, ζ, ι, R, aΔ, bΔ)
    n = size(X,1)
    τ²_n = update_τ²(X, y, μ, γ, Λ, u, D, V)
    u_n, ξ_n = update_u_ξ(u, γ, D, τ²_n, Δ, M, Λ, V)
    γ_n = update_γ(X, D, Λ, u, μ, τ²_n, n)
    D_n = update_D(γ_n, u_n, Λ, θ, τ²_n, V)
    θ_n = update_θ(ζ, ι, V, D_n)
    Δ_n = update_Δ(aΔ, bΔ, ξ, V)
    M_n = update_M(u_n, Λ, V)
    μ_n = update_μ(X, y, γ_n, τ²_n, n)
    Λ_n = update_Λ(πᵥ, R, Λ, u_n, D_n, τ²_n, γ)
    πᵥ_n = update_π(Λ_n, η)
    return (τ²_n, u_n, ξ_n, γ_n, D_n, θ_n, Δ_n, M_n, μ_n, Λ_n, πᵥ_n)
end


function BayesNet(X::Array, y::Array, R::Real; η::Real=1.01,ζ::Real=1,ι::Real=1,aΔ::Real=0,bΔ::Real=0, nburn::Int64=30000, nsamples::Int64=20000, V_in::Int64=NaN, x_transform::Bool=true)
    X, θ, D, πᵥ, Λ, Δ, ξ, M, u, μ, τ², γ, V = init_vars(X, η, ζ, ι, R, aΔ, bΔ, V_in, x_transform)
    
    # burn-in
    for i in 1:nburn
        #τ²[1], u[1], ξ[1], γ[1], D[1], θ[1], Δ[1], M[1], μ[1], Λ[1], πᵥ[1] = GibbsSample(X, y, last(θ), last(D), last(πᵥ), last(Λ), last(Δ), last(ξ), last(M), last(u), last(μ), last(γ), V, η, ζ, ι, R, aΔ, bΔ)
        #TODO: better solution than just flattening last(γ)
        result = GibbsSample(X, y, last(θ), last(D), last(πᵥ), last(Λ), last(Δ), last(ξ), last(M), last(u), last(μ), vec(last(γ)), V, η, ζ, ι, R, aΔ, bΔ)
        push!(τ²,result[1])
        push!(u,result[2])
        push!(ξ,result[3])
        γ = vcat(γ,[result[4]])
        push!(D,result[5])
        push!(θ,result[6])
        push!(Δ,result[7])
        push!(M,result[8])
        push!(μ,result[9])
        push!(Λ,result[10])
        push!(πᵥ,result[11])
    end
    for i in 1:nsamples 
        result = GibbsSample(X, y, last(θ), last(D), last(πᵥ), last(Λ), last(Δ), last(ξ), last(M), last(u), last(μ), vec(last(γ)), V, η, ζ, ι, R, aΔ, bΔ)
        push!(τ²,result[1])
        push!(u,result[2])
        push!(ξ,result[3])
        push!(γ,result[4])
        push!(D,result[5])
        push!(θ,result[6])
        push!(Δ,result[7])
        push!(M,result[8])
        push!(μ,result[9])
        push!(Λ,result[10])
        push!(πᵥ,result[11])
    end
    return τ², u, ξ, γ, D, θ, Δ, M, μ, Λ, πᵥ
end

#main()

