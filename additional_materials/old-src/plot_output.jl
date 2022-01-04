using Plots,Statistics
include("utils.jl")


function plot_γs_test(γ₁, γ₂, γ₁_label, γ₂_label,save_append)
    len = size(γ₁,1)
    lower_bound = convert(Int,floor(len * 0.05))
    upper_bound = convert(Int,ceil(len * 0.95))
    γ₁_sorted = sort(γ₁[:,1])
    γ₂_sorted = sort(γ₂[:,1])
    intervals = scatter([γ₁_sorted[lower_bound],median(γ₁_sorted),γ₁_sorted[upper_bound]],[1,1,1],mark=(:vline),color="red",legend=false, size=(1200,1600))
    plot!(intervals,[γ₁_sorted[lower_bound],γ₁_sorted[upper_bound]],[1,1],color="red",legend=false)
    plot!(intervals,[γ₂_sorted[lower_bound],γ₂_sorted[upper_bound]],[1,1],color="blue",legend=false)
    for i in 2:size(γ₁,2)
        γ₁_sorted = sort(γ₁[:,i])
        γ₂_sorted = sort(γ₂[:,i])
        scatter!(intervals,[γ₁_sorted[lower_bound],median(γ₁_sorted),γ₁_sorted[upper_bound]],[i,i,i],mark=(:vline),color="red",legend=false)
        plot!(intervals,[γ₁_sorted[lower_bound],γ₁_sorted[upper_bound]],[i,i],color="red",legend=false)
        plot!(intervals,[γ₂_sorted[lower_bound],γ₂_sorted[upper_bound]],[i,i],color="blue",legend=false)
    end
    #legend((γ₁_label,γ₂_label))
    savefig(intervals,"plots/test/gamma_cis_$save_append")
end


function plot_γ_sim(γ,title,save_append,jcon)
    len = size(γ,1)
    lower_bound = convert(Int,floor(len * 0.05))
    if lower_bound == 0
        lower_bound = 1
    end
    upper_bound = convert(Int,ceil(len * 0.95))
    γ_sorted = sort(γ[:,1])
    intervals = scatter([γ_sorted[lower_bound],median(γ_sorted),γ_sorted[upper_bound]],[1,1,1],mark=(:vline),color="red",legend=false, size=(1200,1600))
    plot!(intervals,[γ_sorted[lower_bound],γ_sorted[upper_bound]],[1,1],color="red",legend=false)
    for i in 2:size(γ,2)
        γ_sorted = sort(γ[:,i])
        scatter!(intervals,[γ_sorted[lower_bound],median(γ_sorted),γ_sorted[upper_bound]],[i,i,i],mark=(:vline),color="red",legend=false)
        plot!(intervals,[γ_sorted[lower_bound],γ_sorted[upper_bound]],[i,i],color="red",legend=false)
    end
    #legend((γ₁_label,γ₂_label))
    if jcon
        savefig(intervals,"juliacon/plots/gamma_cis_$save_append")
    else
        savefig(intervals,"plots/simulation/gamma_cis_$save_append")
    end
end

function plot_sig_γ(γ,γ₀,simnum,casenum)
    q = size(γ,1)
    V = convert(Int,(1 + sqrt(1 + 8*q))/2)

    selected_idx = select_edges(γ,0.05)
    #selected_idx_true = select_edges(γ₀,0.05)

    sel_γ = zeros(size(γ,1))
    #sel_γ_true = zeros(size(γ₀,1))

    sel_γ[selected_idx] .= 1
    #sel_γ_true[selected_idx_true] .= 1



    sel_B = create_upper_tri(sel_γ,V)
    #sel_B_true = create_upper_tri(sel_γ_true,V)
    plt = plot(sel_B,st=:heatmap,c=palette([:red,:white],10),yflip=true,colorbar=false,
               legend=false,xaxis=(1:V),yaxis=(1:V),size=(1800,1800))

    coords_x = zeros(V,V)
    coords_y = zeros(V,V)
    for i in 1:V
        coords_x[i,:] = 1:V
        coords_y[:,i] = 1:V
    end

    scatter!(plt,upper_triangle(coords_x),upper_triangle(coords_y),series_annotations=round.(γ₀,digits=2),marker=(0,0,0))
    hline!(plt,0.5:(V+0.5),c=:black)
    vline!(plt,0.5:(V+0.5),c=:black)

    savefig(plt,"juliacon/plots/simulation$(simnum)_case$(casenum)_posterior_B")
    #savefig(plot(sel_B_true,st=:heatmap),"juliacon/plots/simulation$(simnum)_case$(casenum)_true_B")

end

function plot_ξ_sim(ξ_df, simnum, casenum)
    plt = bar(ξ_df[:,"Xi posterior"],fillcolor=ξ_df[:,"Color"],legend=false,xlabel="Microbe",ylabel="Probability of influence")

    savefig(plt,"juliacon/plots/simulation$(simnum)_case$(casenum)_posterior_xi")
end
