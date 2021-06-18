using Plots

function plot_γs_test(γ₁, γ₂, γ₁_label, γ₂_label,save_append)
    len = size(γ₁,2)
    lower_bound = convert(Int,floor(len * 0.05))
    upper_bound = convert(Int,ceil(len * 0.95))
    γ₁_sorted = sort(γ₁[1,:])
    γ₂_sorted = sort(γ₂[:,1])
    intervals = scatter([γ₁_sorted[lower_bound],median(γ₁_sorted),γ₁_sorted[upper_bound]],[1,1,1],mark=(:vline),color="red",legend=false, size=(1200,1600))
    plot!(intervals,[γ₁_sorted[lower_bound],γ₁_sorted[upper_bound]],[1,1],color="red",legend=false)
    plot!(intervals,[γ₂_sorted[lower_bound],γ₂_sorted[upper_bound]],[1,1],color="blue",legend=false)
    for i in 2:190
        γ₁_sorted = sort(γ₁[i,:])
        γ₂_sorted = sort(γ₂[:,i])
        scatter!(intervals,[γ₁_sorted[lower_bound],median(γ₁_sorted),γ₁_sorted[upper_bound]],[i,i,i],mark=(:vline),color="red",legend=false)
        plot!(intervals,[γ₁_sorted[lower_bound],γ₁_sorted[upper_bound]],[i,i],color="red",legend=false)
        plot!(intervals,[γ₂_sorted[lower_bound],γ₂_sorted[upper_bound]],[i,i],color="blue",legend=false)
    end
    #legend((γ₁_label,γ₂_label))
    savefig(intervals,"plots/test/gamma_cis_$save_append")
end


function plot_γ_sim(γ,title,save_append)
    len = size(γ,2)
    lower_bound = convert(Int,floor(len * 0.05))
    upper_bound = convert(Int,ceil(len * 0.95))
    γ_sorted = sort(γ[1,:])
    intervals = scatter([γ_sorted[lower_bound],median(γ_sorted),γ_sorted[upper_bound]],[1,1,1],mark=(:vline),color="red",legend=false, size=(1200,1600))
    plot!(intervals,[γ_sorted[lower_bound],γ_sorted[upper_bound]],[1,1],color="red",legend=false)
    for i in 2:190
        γ_sorted = sort(γ[i,:])
        scatter!(intervals,[γ_sorted[lower_bound],median(γ_sorted),γ_sorted[upper_bound]],[i,i,i],mark=(:vline),color="red",legend=false)
        plot!(intervals,[γ_sorted[lower_bound],γ_sorted[upper_bound]],[i,i],color="red",legend=false)
    end
    #legend((γ₁_label,γ₂_label))
    savefig(intervals,"plots/simulation/gamma_cis_$save_append")
end
