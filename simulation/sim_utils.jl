using DrWatson,DataFrames,CSV

function fill_B(B,ξ,t,μₛ,σₛ,rng)
    for i in 1:t
        for j in (i+1):t
            if (ξ[i] == 1 && ξ[j] == 1)
                B[i,j] = B[j,i] = rand(rng,Normal(μₛ,σₛ))
            end
        end
    end
end

function output_data(saveinfo,out_df,B,m,ξ,jcon,simtype)

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

end