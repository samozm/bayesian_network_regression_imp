using Random, Distributions, DataFrames, LinearAlgebra, StatsBase, RCall, InvertedIndices
R"library(GIGrvg)"

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
        wts[0] = r^η + 1
    elseif(λ[r] == 0)
        wts[1] = 2
    else
        wts[2] = 2
    end
    rand(Dirichlet(wts))
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
    - `S` : set to a matrix of V×V draws from the Exponential(θ/2) distribution
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
    X_new = transpose(hcat(map(k -> upper_triangle(X[k]), 1:size(X,1))...))

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
    ξ = map(keys -> rand(Binomial(1,Δ)), 1:V)
    M = rand(InverseWishart(V,Matrix(I,R,R)))
    u = hcat(map(k -> sample_u(ξ[k], R, M), 1:V)...)
    μ = 1
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
    yμ1Xγ = (y - μ*ones(n,1) - X*γ)
    γW = (γ - W)
    yμ1Xγᵀyμ1Xγ = transpose(yμ1Xγ)*yμ1Xγ
    γWᵀγW = transpose(γW)*inv(D)*γW

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
    D = diagm(map(k -> sample_rgig(a[k],θ), 1:q))
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
        U = u[:,Not(k)]*Λ
        s = create_upper_tri(diag(D),V)
        H = diagm(vcat(filter!(i->i!=0,s[:,k]),filter!(i->i!=0,s[k,:])))
        Γ = create_upper_tri(γ, V)
        γk= vcat(filter!(i->i!=0,Γ[:,k]),filter!(i->i!=0,Γ[k,:]))
        Σ = inv(((transpose(U)*inv(H)*U)/τ²) + inv(M))
        m = (Σ*transpose(U)*inv(H)*γk)/τ²
        mvn_a = MultivariateNormal(zeros(size(H,1)),τ²*H)
        mvn_b_Σ = round.(τ²*H + U*M*transpose(U), digits=10)
        mvn_b = MultivariateNormal(zeros(size(H,1)),mvn_b_Σ)
        w_top = (1-Δ) * pdf(mvn_a,γk)
        w_bot = w_top + Δ*pdf(mvn_b,γk)
        w = w_top / w_bot
        mvn_f = MultivariateNormal(m,round.(Σ,digits=10))
        
        ξ[k] = update_ξ(w)
        # the paper says take the pdf of mvn_f at u_k, but that doesn't really make sense since we need R values not 1
        # in their implementation they sample from mvn_f
        # also, the paper says the first term is (1-w) but their code uses 1-ξ. again i think this makes more sense
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
    b = bΔ + V - sum(ξ)
    return rand(Beta(a,b))
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
        # in the paper this is uᵀΛu but in the code it's just uuᵀ
        uΛu = uΛu .+ u[:,i]*transpose(u[:,i])#*Λ*u[:,i]
        if u[:,i] != zeros(size(u,1))
            num_nonzero = num_nonzero + 1
        end
    end
    println("uΛu")
    show(stdout, "text/plain", uΛu)
    println("")
    println("I(Returns)")
    show(stdout, "text/plain", I(R))
    println("")
    Ψ = I(R) .+ uΛu
    df = V + num_nonzero
    println("Ψ ")
    show(stdout, "text/plain", Ψ )
    println("")
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

# Returns
new value of Λ
"""
function update_Λ(πᵥ, R, Λ, u, D, τ²)
    λ = diag(Λ)
    λ_new = zeros(size(Λ,1))
    for r in 1:R
        Λ₋₁= diagm(hcat(λ[1:r-1],[-1],λ[r+1:R]))
        println("Λ₋₁")
        show(stdout, "text/plain", Λ₋₁)
        println("")
        Λ₀ = diagm(hcat(λ[1:r-1],[0],λ[r+1:R]))
        Λ₁ = diagm(hcat(λ[1:r-1],[1],λ[r+1:R]))
        W₋₁= upper_triangle(u * Λ₋₁ * transpose(u))
        W₀ = upper_triangle(u * Λ₀ * transpose(u))
        W₁ = upper_triangle(u * Λ₁ * transpose(u))
        n₀ = pdf(MultivariateNormal(transpose(W₀), τ² * D),γ)
        n₁ = pdf(MultivariateNormal(transpose(W₁), τ² * D),γ)
        n₋₁= pdf(MultivariateNormal(transpose(W₋₁), τ² * D),γ)
        p_bot = πᵥ[:,1] * n₀ + πᵥ[:,2] * n₁ + πᵥ[:3] * n₋₁
        p1 = πᵥ[:,1] * n₀ / p_bot
        p2 = πᵥ[:,2] * n₁ / p_bot
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
    Λ_n = update_Λ(πᵥ, R, Λ, u_n, D_n, τ²_n)
    πᵥ_n = update_π(Λ_n, η)
    return (τ²_n, u_n, ξ_n, γ_n, D_n, θ_n, Δ_n, M_n, μ_n, Λ_n, πᵥ_n)
end


function BayesNet(X::Array{Array{T,2},1}, y::Array, R::Real,η::Real=1.01,ζ::Real=1,ι::Real=1,aΔ::Real=0,bΔ::Real=0, nburn::Int64=30000, nsamples::Int64=20000) where T <: Real
    X, θ, D, πᵥ, Λ, Δ, ξ, M, u, μ, τ², γ, V = init_vars(X, η, ζ, ι, R, aΔ, bΔ)
    
    # burn-in
    for i in 1:nburn
        τ²[1], u[1], ξ[1], γ[1], D[1], θ[1], Δ[1], M[1], μ[1], Λ[1], πᵥ[1] = GibbsSample(X, y, last(θ), last(D), last(πᵥ), last(Λ), last(Δ), last(ξ), last(M), last(u), last(μ), last(γ), V, η, ζ, ι, R, aΔ, bΔ)
    end
    for i in 1:nsamples 
        result = GibbsSample(X, y, last(θ), last(D), last(πᵥ), last(Λ), last(Δ), last(ξ), last(M), last(u), last(μ), last(γ), V, η, ζ, ι, R, aΔ, bΔ)
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

