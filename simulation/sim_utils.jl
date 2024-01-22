using DrWatson,DataFrames,CSV,Statistics

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


function augment(X,y,rng)
    keep_prob = 0.9
    new_sd = Statistics.std(y)

    rbin = Bernoulli(keep_prob)
    rnorm = Normal(0,new_sd)
    rchisq = Chisq(3)/15

    X_new = zeros(size(X,1)*2,size(X,2))
    y_new = zeros(size(y,1)*2)
    LEN_OTU = size(X,1)
    
    for i in 1:LEN_OTU
        num_necessary = size(X,2)
        for j in 1:num_necessary
            X_new[i,j] = X[i,j]
            X_new[i+LEN_OTU,j] = ifelse(X[i,j] > 0, rand(rng,rbin), rand(rng,rbin))
        end

        y_new[i] = y[i]
        y_new[i+LEN_OTU] = y[i] + rand(rng,rnorm)
        if y_new[i+LEN_OTU] <= 0
            y_new[i+LEN_OTU] = rand(rng,rchisq)
        end
    end
    return X_new,y_new
end