using GaussianMixtures

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
