using DataFrames,CSV

unreal_psrf = DataFrame()
real_psrf = DataFrame()
normed_psrf = DataFrame()

for flnm in readdir("results/simulation/local/unrealistic-results")
    if !occursin("psrf",flnm)
        continue
    end
    mu = flnm[last(findfirst("mu=",flnm))+1:last(findfirst("mu=",flnm))+3]
    pi_ = flnm[last(findfirst("pi=",flnm))+1:last(findfirst("pi=",flnm))+3]
    R = flnm[last(findfirst("R=",flnm))+1]
    k = flnm[last(findfirst("n_microbes=",flnm))+1]
    if k == '2'
        k = 22
    elseif k == '1'
        k = 15
    elseif k == '8'
        k = 8
    elseif k == '3'
        k = 30
    else
        print("bad switch")
    end
    n = flnm[last(findfirst("samplesize=",flnm))+1:last(findfirst("samplesize=",flnm))+3]

    df_tmp = DataFrame(CSV.File("results/simulation/local/unrealistic-results/$flnm"))

    append!(unreal_psrf,DataFrame(mu=mu,pi=pi_,R=R,k=k,n=n,xi=round(df_tmp[1,:max_xi],digits=2),gamma=round(df_tmp[1,:max_gamma],digits=2)))
end

for flnm in readdir("results/simulation/local/realistic-results")
    if !occursin("out=psrf",flnm)
        continue
    end
    if occursin("psrf=1.",flnm)
        continue
    end
    mu = flnm[last(findnext("mu=",flnm,14))+1:last(findnext("mu=",flnm,14))+3]
    pi_ = flnm[last(findfirst("pi=",flnm))+1:last(findfirst("pi=",flnm))+3]
    R = flnm[last(findfirst("R=",flnm))+1]
    k = flnm[last(findfirst("n_microbes=",flnm))+1]
    nu = flnm[last(findfirst("nu=",flnm))+1:first(findfirst("_out=",flnm))-1]
    if k == '2'
        k = 22
    elseif k == '1'
        k = 15
    elseif k == '8'
        k = 8
    elseif k == '3'
        k = 30
    else
        print("bad switch")
    end
    if R <= '3'
        continue
    end
    n = flnm[last(findfirst("samplesize=",flnm))+1:first(findfirst("_simnum=",flnm))-1]
    if (n == "50") && (nu != "12")
        continue
    end

    ix1 = last(findfirst("type=",flnm))+1
    ix2 = first(findfirst(".csv",flnm))-1
    type=flnm[ix1:ix2]
    out=type

    df_tmp = DataFrame(CSV.File("results/simulation/local/realistic-results/$flnm"))

    append!(real_psrf,DataFrame(mu=mu,pi=pi_,R=R,k=k,n=n,nu=nu,type=type,xi=round(df_tmp[1,:max_xi],digits=2),gamma=round(df_tmp[1,:max_gamma],digits=2)))
end

for flnm in readdir("results/simulation/local/normed-results")
    if !occursin("psrf",flnm)
        continue
    end
    mu = flnm[last(findnext("mu=",flnm,14))+1:last(findnext("mu=",flnm,14))+3]
    pi_ = flnm[last(findfirst("pi=",flnm))+1:last(findfirst("pi=",flnm))+3]
    R = flnm[last(findfirst("R=",flnm))+1]
    k = flnm[last(findfirst("n_microbes=",flnm))+1]
    if k == '2'
        k = 22
    elseif k == '1'
        k = 15
    elseif k == '8'
        k = 8
    elseif k == '3'
        k = 30
    else
        print("bad switch")
    end
    n = flnm[last(findfirst("samplesize=",flnm))+1:first(findfirst("_simnum=",flnm))-1]

    ix1 = last(findfirst("type=",flnm))+1
    ix2 = first(findfirst(".csv",flnm))-1
    type=flnm[ix1:ix2]
    out=type

    df_tmp = DataFrame(CSV.File("results/simulation/local/normed-results/$flnm"))

    append!(normed_psrf,DataFrame(mu=mu,pi=pi_,R=R,k=k,n=n,type=type,xi=df_tmp[1,:max_xi],gamma=df_tmp[1,:max_gamma]))
end

println("")
println("Unrealistic > 1.01")
println("")

display(filter(row -> (row.gamma > 1.01 || row.xi > 1.01), unreal_psrf))

println("")
println("")
println("--------- Realistic ----------")


println("")
println("Realistic > 1.01")
println("")

display(filter(row -> (row.gamma > 1.01 || row.xi > 1.01), real_psrf)) 

println("")
println("Realistic > 1.03")
println("")

display(filter(row -> (row.gamma > 1.03 || row.xi > 1.03), real_psrf)) 

println("")
println("Realistic > 1.1")
println("")

display(filter(row -> (row.gamma > 1.1 || row.xi > 1.1), real_psrf)) 
println("")
println("")

#println("")
#println("--------- Normed ----------")


#println("")
#println("Normed > 1.01")
#println("")

#display(filter(row -> (row.gamma > 1.01 || row.xi > 1.01), normed_psrf)) 

#println("")
#println("Normed > 1.05")
#println("")

#display(filter(row -> (row.gamma > 1.05 || row.xi > 1.05), normed_psrf)) 

#println("")
#println("Normed > 1.1")
#println("")

#display(filter(row -> (row.gamma > 1.1 || row.xi > 1.1), normed_psrf)) 

#println("")


edges = DataFrame()

for flnm in readdir("results/simulation/local/realistic-results")
    if !occursin("edges",flnm)
        continue
    end
    if !occursin("psrf=1.1",flnm)
        continue
    end
    mu = flnm[last(findnext("mu=",flnm,14))+1:last(findnext("mu=",flnm,14))+3]
    pi_ = flnm[last(findfirst("pi=",flnm))+1:last(findfirst("pi=",flnm))+3]
    R = flnm[last(findfirst("R=",flnm))+1]
    k = flnm[last(findfirst("n_microbes=",flnm))+1]
    if k == '2'
        k = 22
    elseif k == '1'
        k = 15
    elseif k == '8'
        k = 8
    elseif k == '3'
        k = 30
    else
        print("bad switch")
    end
    n = flnm[last(findfirst("samplesize=",flnm))+1:first(findfirst("_simnum=",flnm))-1]

    ix1 = last(findfirst("type=",flnm))+1
    ix2 = first(findfirst(".csv",flnm))-1
    type=flnm[ix1:ix2]
    out=type

    df_tmp = DataFrame(CSV.File("results/simulation/local/realistic-results/$flnm"))

    append!(edges,DataFrame(mu=mu,pi=pi_,R=R,k=k,n=n,type=type,mean=df_tmp[1,:mean]))
end

#println("")
#println("--------- Edges ----------")


#println("")
#println("|Edges| > 7")
#println("")

#display(filter(row -> (row.mean > 7 || row.mean < -7), edges))
#println("")