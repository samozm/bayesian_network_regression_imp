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

function output_data(saveinfo,out_df,ξ,A,B,jcon,simtype)

    if jcon
        saveinfo["out"] = "XYs"
        CSV.write(projectdir("juliacon","data",savename(saveinfo,"csv",digits=1)),out_df)

        saveinfo["out"] = "xis"
        CSV.write(projectdir("juliacon","data",savename(saveinfo,"csv",digits=1)),DataFrame(TrueXi=ξ))

        saveinfo["out"] = "bs"
        CSV.write(datadir(joinpath("simulation",simtype),savename(saveinfo,"csv",digits=1)),DataFrame(B=lower_triangle(B)[:,1]))
        return
    end
    saveinfo["out"] = "XYs"
    CSV.write(datadir(joinpath("simulation",simtype),savename(saveinfo,"csv",digits=1)),out_df)
    
    saveinfo["out"] = "xis"
    CSV.write(datadir(joinpath("simulation",simtype),savename(saveinfo,"csv",digits=1)),DataFrame(TrueXi=ξ))

    saveinfo["out"] = "A"
    CSV.write(datadir(joinpath("simulation",simtype),savename(saveinfo,"csv",digits=1)),DataFrame(A,:auto))

    saveinfo["out"] = "bs"
    CSV.write(datadir(joinpath("simulation",simtype),savename(saveinfo,"csv",digits=1)),B)

end

function augment(X,y,ξ,A_base,rng,num_aug)
    keep_prob = 0.9
    new_sd = Statistics.std(y)

    rbin = Bernoulli(keep_prob)
    rnorm = Normal(0,new_sd)
    rchisq = Chisq(3)/15

    num_aug = convert(Int64,num_aug)

    A_X = lower_triangle(A_base)

    X_new = zeros(size(X,1)*(num_aug+1),size(X,2))
    y_new = zeros(size(y,1)*(num_aug+1))
    LEN_OTU = size(X,1)
    LEN_XI = size(ξ,1)
    
    for xx in 1:num_aug
        for i in 1:LEN_OTU
            num_necessary = size(X,2)
            for j in 1:num_necessary
                X_new[i,j] = X[i,j]
                X_new[i+(LEN_OTU*xx),j] = ifelse(X_new[i,j] > 0, rand(rng,rbin) * A_X[j], (1-rand(rng,rbin)) * A_X[j])
            end

            y_new[i] = y[i]
            y_new[i+(LEN_OTU*xx)] = y[i] + rand(rng,rnorm)
            #if y_new[i+(LEN_OTU*xx)] <= 0
            #    y_new[i+(LEN_OTU*xx)] = rand(rng,rchisq)
            #end
        end
    end
    return X_new,y_new
end

# from https://atrebas.github.io/post/2021-01-17-index_to_lower_triangular_subscripts/
function lin_idx_to_lt_idx(i::Int,sz)
    p = (sqrt(1 + 8 * i) - 1)/2
    i₀ = Int64(floor(p))
    if i₀ == p
        return i₀,i₀
    else
        return i₀ + 1, Int64(i - i₀*(i₀+1)/2)
    end
end


function create_lower_tri(vector::AbstractVector{T},V) where {T}
    mat = zeros(T,V,V)
    i = 1
    for k = 1:V
        for l = k:V
            mat[l,k] = vector[i,1] 
            i += 1
        end
    end
    return mat
end

function lower_triangle(matrix::AbstractArray{T,2}) where {T}
    if size(matrix,1) != size(matrix,2)
        print("error: matrix must be square")
        return 0
    end

    V = size(matrix,1)

    k = 1
    ret = zeros(T,convert(Int64, round(V*(V + 1)/2)))
    for i in 1:size(matrix,1)
        for j in i:size(matrix,2)
            ret[k,1] = matrix[j,i]
            k = k + 1
        end
    end
    return ret
end