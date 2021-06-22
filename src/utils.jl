using GaussianMixtures


### implementing FDR controlling edge selection as described in appendix c
function select_edges(means,α)
    abslog_means = DataFrame(idx=1:size(means,1),means=log.(abs.(means))
    gmm = GMM(2,abslog_means[:means])
    pb = gmmposterior(gmm)[1]
    gmm_means = means(gmm)
    fdr = 0
    H = 1
    sorted_abslog_means = sort(abslog_means,:means,rev=true)

    #TODO: find largest H with FDR lower than alpha
    while true
        if H > 1
            new_fdr = (fdr*(H-1) + pb[sorted_abslog_means[H,:idx],argmin(gmm_means)])/H
        else
            new_fdr = pb[sorted_abslog_means[H,:idx],argmin(gmm_means)]
        if new_fdr > α
            break
        end
        H = H+1
        fdr = new_fdr
    end
    return(sorted_abslog_means[1:H,:idx])
end
