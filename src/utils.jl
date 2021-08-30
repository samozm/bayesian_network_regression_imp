using GaussianMixtures,BayesianNetworkRegression,ProgressMeter,MCMCChains,TypedTables
using Random,Distributions

"""
    create_upper_tri(vec,V)

Create an upper triangluar matrix from a vector of the form [12, ... 1V,23,...(V-1)V]
to the form [0 12 13  ... 1V]
            [0 0  23  ... 2V]
            [...............]
            [0 0  0...(V-1)V]
# Arguments
- `vec`: vector containing values to put into the upper triangluar matrix
- `V`  : dimension of output matrix

# Returns
Upper triangluar matrix containing values of `vec`
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

"""
    upper_triangle(matrix)

Return the upper triangle (without the diagonal) of the matrix as a vector

# Arguments
- `matrix`: matrix of which to capture the upper triangle

# Returns
Vector of upper triangluar section of `matrix`
"""
function upper_triangle(matrix)
    k = 1
    ret = zeros(convert(Int64, round(size(matrix,1)*(size(matrix,2) - 1)/2)))
    for i in 1:size(matrix,1)
        for j in (i+1):size(matrix,2)
            ret[k] = matrix[i,j]
            k = k + 1
        end
    end
    return ret
end

function gelman_rubin(X::AbstractArray{T,2}, y::AbstractVector{U}, nburn, nsamples) where {T,U}

    R = 9
    n = size(X,1)
    q = size(X,2)
    total = nburn + nsamples + 1
    V = convert(Int,(1 + sqrt(1 + 8*q))/2)

    nstate = 3
    
    state1 = Table(τ² = Array{Float64,3}(undef,(total,1,1)), u = Array{Float64,3}(undef,(total,R,V)),
                ξ = Array{Float64,3}(undef,(total,V,1)), γ = Array{Float64,3}(undef,(total,q,1)),
                S = Array{Float64,3}(undef,(total,q,1)), θ = Array{Float64,3}(undef,(total,1,1)),
                Δ = Array{Float64,3}(undef,(total,1,1)), M = Array{Float64,3}(undef,(total,R,R)),
                μ = Array{Float64,3}(undef,(total,1,1)), λ = Array{Float64,3}(undef,(total,R,1)),
                πᵥ= Array{Float64,3}(undef,(total,R,3)))
    state2 = Table(τ² = Array{Float64,3}(undef,(total,1,1)), u = Array{Float64,3}(undef,(total,R,V)),
                ξ = Array{Float64,3}(undef,(total,V,1)), γ = Array{Float64,3}(undef,(total,q,1)),
                S = Array{Float64,3}(undef,(total,q,1)), θ = Array{Float64,3}(undef,(total,1,1)),
                Δ = Array{Float64,3}(undef,(total,1,1)), M = Array{Float64,3}(undef,(total,R,R)),
                μ = Array{Float64,3}(undef,(total,1,1)), λ = Array{Float64,3}(undef,(total,R,1)),
                πᵥ= Array{Float64,3}(undef,(total,R,3)))
    state3 = Table(τ² = Array{Float64,3}(undef,(total,1,1)), u = Array{Float64,3}(undef,(total,R,V)),
                ξ = Array{Float64,3}(undef,(total,V,1)), γ = Array{Float64,3}(undef,(total,q,1)),
                S = Array{Float64,3}(undef,(total,q,1)), θ = Array{Float64,3}(undef,(total,1,1)),
                Δ = Array{Float64,3}(undef,(total,1,1)), M = Array{Float64,3}(undef,(total,R,R)),
                μ = Array{Float64,3}(undef,(total,1,1)), λ = Array{Float64,3}(undef,(total,R,1)),
                πᵥ= Array{Float64,3}(undef,(total,R,3)))
    #state4 = Table(τ² = Array{Float64,3}(undef,(total,1,1)), u = Array{Float64,3}(undef,(total,R,V)),
    #            ξ = Array{Float64,3}(undef,(total,V,1)), γ = Array{Float64,3}(undef,(total,q,1)),
    #            S = Array{Float64,3}(undef,(total,q,1)), θ = Array{Float64,3}(undef,(total,1,1)),
    #            Δ = Array{Float64,3}(undef,(total,1,1)), M = Array{Float64,3}(undef,(total,R,R)),
    #            μ = Array{Float64,3}(undef,(total,1,1)), λ = Array{Float64,3}(undef,(total,R,1)),
    #            πᵥ= Array{Float64,3}(undef,(total,R,3)))
    #state5 = Table(τ² = Array{Float64,3}(undef,(total,1,1)), u = Array{Float64,3}(undef,(total,R,V)),
    #            ξ = Array{Float64,3}(undef,(total,V,1)), γ = Array{Float64,3}(undef,(total,q,1)),
    #            S = Array{Float64,3}(undef,(total,q,1)), θ = Array{Float64,3}(undef,(total,1,1)),
    #            Δ = Array{Float64,3}(undef,(total,1,1)), M = Array{Float64,3}(undef,(total,R,R)),
    #            μ = Array{Float64,3}(undef,(total,1,1)), λ = Array{Float64,3}(undef,(total,R,1)),
    #            πᵥ= Array{Float64,3}(undef,(total,R,3)))

    X_new = Matrix{eltype(T)}(undef, n, q)
    #BayesianNetworkRegression.initialize_variables!(state1,X_new,X,1.01,1,0.1,9,1,1,10,V,false)
    #BayesianNetworkRegression.initialize_variables!(state2,X_new,X,1.01,1,0.1,9,1,12,14,V,false)
    #BayesianNetworkRegression.initialize_variables!(state3,X_new,X,1.01,1,0.6,9,8,1,20,V,false)
    #BayesianNetworkRegression.initialize_variables!(state4,X_new,X,1.01,1,0.6,9,8,12,26,V,false)

    BayesianNetworkRegression.initialize_variables!(state1,X_new,X,1.01,1,1,9,1,1,10,V,false)
    BayesianNetworkRegression.initialize_variables!(state2,X_new,X,1.01,1,1,9,1,1,10,V,false)
    BayesianNetworkRegression.initialize_variables!(state3,X_new,X,1.01,1,1,9,1,1,10,V,false)
    #BayesianNetworkRegression.initialize_variables!(state4,X_new,X,1.01,1,1,9,1,1,10,V,false)
    #BayesianNetworkRegression.initialize_variables!(state5,X_new,X,1.01,1,1,9,1,1,10,V,false)
    for i in 1:q
        state1.γ[1,i,1] = rand(Normal(0,8))
        state2.γ[1,i,1] = rand(Normal(0,8))
        state3.γ[1,i,1] = rand(Normal(0,8))
        #state4.γ[1,i,1] = rand(Normal(0,3))
        #tate5.γ[1,i,1] = rand(Normal(0,3))
    end
    
    state1.ξ[1,:,1] = vcat(repeat([0],Int64(V)))
    state2.ξ[1,:,1] = vcat(repeat([1],Int64(V)))
    state3.ξ[1,:,1] = repeat([0,1],Int64(V/2))
    #state4.ξ[1,:,1] = repeat([1,0,0,0,1,0,0,0,1,0,0,0,1,0,0],2)
    #state5.ξ[1,:,1] = repeat([0,1,1,1,0,1,1,1,0,1,1,1,0,1,1],2)

    p = Progress(total-1,1)
    for i in 2:total
        BayesianNetworkRegression.GibbsSample!(state1, i, X_new, y, V, 1.01, 1, 1, 9, 1, 1, 10)
        BayesianNetworkRegression.GibbsSample!(state2, i, X_new, y, V, 1.01, 1, 1, 9, 1, 1, 10)
        BayesianNetworkRegression.GibbsSample!(state3, i, X_new, y, V, 1.01, 1, 1, 9, 1, 1, 10)
        #BayesianNetworkRegression.GibbsSample!(state4, i, X_new, y, V, 1.01, 1, 1, 9, 1, 1, 10)
        #BayesianNetworkRegression.GibbsSample!(state5, i, X_new, y, V, 1.01, 1, 1, 9, 1, 1, 10)
        next!(p)
    end

    all_γs = Array{Float64,3}(undef,(nsamples,q,nstate))
    all_γs[:,:,1] = state1.γ[nburn+2:total,:,1]
    all_γs[:,:,2] = state2.γ[nburn+2:total,:,1]
    all_γs[:,:,3] = state3.γ[nburn+2:total,:,1]
    #all_γs[:,:,4] = state4.γ[:,:,1]
    #all_γs[:,:,5] = state5.γ[:,:,1]

    all_ξs = Array{Float64,3}(undef,(nsamples,V,nstate))
    all_ξs[:,:,1] = state1.ξ[nburn+2:total,:,1]
    all_ξs[:,:,2] = state2.ξ[nburn+2:total,:,1]
    all_ξs[:,:,3] = state3.ξ[nburn+2:total,:,1]
    #all_ξs[:,:,4] = state4.ξ[:,:,1]
    #all_ξs[:,:,5] = state5.ξ[:,:,1]
    
    #return MCMCChains.gelmandiag_multivariate(all_γs)
    return MCMCChains.gelmandiag(all_γs),MCMCChains.gelmandiag(all_ξs)

end

function gelman_rubin_intialize_values(τ²,u,ξ,γ,S,θ,Δ,M,μ,λ,πᵥ,R,total,V,q)
    state = Table(τ² = Array{Float64,3}(undef,(total,1,1)), u = Array{Float64,3}(undef,(total,R,V)),
                ξ = Array{Float64,3}(undef,(total,V,1)), γ = Array{Float64,3}(undef,(total,q,1)),
                S = Array{Float64,3}(undef,(total,q,1)), θ = Array{Float64,3}(undef,(total,1,1)),
                Δ = Array{Float64,3}(undef,(total,1,1)), M = Array{Float64,3}(undef,(total,R,R)),
                μ = Array{Float64,3}(undef,(total,1,1)), λ = Array{Float64,3}(undef,(total,R,1)),
                πᵥ= Array{Float64,3}(undef,(total,R,3)))

    state.θ[1] = θ

    state.S[1,:] = S
    state.πᵥ[1,:,:] = πᵥ
    
    state.λ[1,:] = λ
    state.Δ[1] = Δ

    state.ξ[1,:] = ξ
    state.M[1,:,:] = M
    
    state.u[1,:,:] = u
    state.μ[1] = μ

    state.τ²[1] = τ²
    state.γ[1,:] = γ
    
    return state
end


### implementing FDR controlling edge selection as described in appendix c
function select_edges(γ_means,α)
    abslog_means = DataFrame(idx=1:size(γ_means,1),means=log.(abs.(γ_means)))
    #abslog_means = DataFrame(idx=1:size(γ_means,1),means=abs.(γ_means))
    #println("")
    #show(stdout,"text/plain",abslog_means[:,:means])
    #println("")
    gmm = GMM(2,reshape(abslog_means[:,:means],(size(γ_means,1),1)))
    #gmm = GMM(2,1,kind=:full)
    #em!(gmm,reshape(abslog_means[:,:means],(size(γ_means,1),1)))
    pb = gmmposterior(gmm,reshape(abslog_means[:,:means],(size(γ_means,1),1)))[1]
    gmm_means = means(gmm)
    fdr = 0
    new_fdr = 0
    H = 1
    sorted_abslog_means = abslog_means
    sorted_abslog_means[:,"origmeans"] = γ_means
    sorted_abslog_means = sort(sorted_abslog_means,:means,rev=true)

    while true
        fdr = new_fdr
        new_fdr = (new_fdr*(H-1) + pb[sorted_abslog_means[H,:idx],argmin(gmm_means)])/H
        #println("fdr")
        #println(pb[sorted_abslog_means[H,:idx],:])
        #println(min(pb[:,argmin(gmm_means)]...))
        #println("")
        #show(stdout,"text/plain",sorted_abslog_means)
        #println("")
        pb2 = DataFrame(pb,:auto)
        pb2[:,"idx"] = 1:size(γ_means,1)
        #show(stdout,"text/plain",pb2)
        #println("")
        #show(stdout,"text/plain",gmm_means)
        #println("H")
        #println(H)
        if new_fdr > α || H == size(γ_means,1)
            break
        end
        H = H+1
    end
    return gmm#sorted_abslog_means[1:H,:idx]
end
