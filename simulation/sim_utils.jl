using DrWatson,DataFrames,CSV

function fill_B(B,ξ,t,μₛ,σₛ)
    for i in 1:t
        for j in (i+1):t
            if (ξ[i] == 1 && ξ[j] == 1)
                B[i,j] = B[j,i] = rand(Normal(μₛ,σₛ))
            end
        end
    end
end

function output_data(saveinfo,out_df,B,m,ξ,A,jcon,simtype)

    if jcon
        saveinfo["out"] = "XYs"
        CSV.write(projectdir("juliacon","data",savename(saveinfo,"csv",digits=1)),out_df)

        saveinfo["out"] = "bs"
        CSV.write(projectdir("juliacon","data",savename(saveinfo,"csv",digits=1)),DataFrame(B=BayesianNetworkRegression.lower_triangle(B)[:,1]))

        saveinfo["out"] = "ms"
        CSV.write(projectdir("juliacon","data",savename(saveinfo,"csv",digits=1)),DataFrame(transpose(hcat(m...)),:auto))

        saveinfo["out"] = "xis"
        CSV.write(projectdir("juliacon","data",savename(saveinfo,"csv",digits=1)),DataFrame(TrueXi=ξ))
        return
    end
    saveinfo["out"] = "XYs"
    CSV.write(datadir(joinpath("simulation",simtype),savename(saveinfo,"csv",digits=1)),out_df)
    
    saveinfo["out"] = "bs"
    CSV.write(datadir(joinpath("simulation",simtype),savename(saveinfo,"csv",digits=1)),DataFrame(B=BayesianNetworkRegression.lower_triangle(B)[:,1]))
   
    saveinfo["out"] = "ms"
    CSV.write(datadir(joinpath("simulation",simtype),savename(saveinfo,"csv",digits=1)),DataFrame(transpose(hcat(m...)),:auto))
    
    saveinfo["out"] = "xis"
    CSV.write(datadir(joinpath("simulation",simtype),savename(saveinfo,"csv",digits=1)),DataFrame(TrueXi=ξ))

    V = size(A[1],1)
    q = Int64(V*(V-1)/2)
    n = size(A,1)
    A_reshaped = Matrix{Float64}(undef, n, q)
    for i in 1:n
        A_reshaped[i,:] = BayesianNetworkRegression.lower_triangle(A[i])
    end

    saveinfo["out"] = "As"
    CSV.write(datadir(joinpath("simulation",simtype),savename(saveinfo,"csv",digits=1)),DataFrame(A_reshaped,:auto))

end