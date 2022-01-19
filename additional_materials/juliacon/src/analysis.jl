using DataFrames,CSV

include("../../old-src/utils.jl")
include("../../old-src/plot_output.jl")


function analyze_γ(simnum,casenum)
    γ = upper_triangle(Matrix(DataFrame(CSV.File("juliacon/results/simulation$(simnum)_case$(casenum)_gammas.csv"))))
    γ₀= DataFrame(CSV.File("juliacon/data/simulation$(simnum)_case$(casenum)_bs.csv"))

    plot_sig_γ(γ,γ₀[:,:B],simnum,casenum)
end

function analyze_ξ(simnum,casenum)
    ξ_df = DataFrame(CSV.File("juliacon/results/simulation$(simnum)_case$(casenum).csv"))
    ξ_df[:,"Color"] .= "Blue"
    ξ_df[ξ_df.TrueXi .== true,"Color"] .= "Red"

    plot_ξ_sim(ξ_df,simnum,casenum)

end
