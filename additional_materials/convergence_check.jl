using DataFrames,CSV

unreal_psrf = DataFrame()
real_psrf = DataFrame()

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
    else
        print("bad switch")
    end
    n = flnm[last(findfirst("samplesize=",flnm))+1:last(findfirst("samplesize=",flnm))+3]

    df_tmp = DataFrame(CSV.File("results/simulation/local/unrealistic-results/$flnm"))

    append!(unreal_psrf,DataFrame(mu=mu,pi=pi_,R=R,k=k,n=n,xi=df_tmp[1,:max_xi],gamma=df_tmp[1,:max_gamma]))
end

for flnm in readdir("results/simulation/local/realistic-results")
    if !occursin("psrf",flnm)
        continue
    end
    mu = flnm[last(findnext("mu=",flnm,14))+1:last(findnext("mu=",flnm,14))+3]
    pi_ = flnm[last(findfirst("pi=",flnm))+1:last(findfirst("pi=",flnm))+3]
    R = flnm[last(findfirst("R=",flnm))+1]
    k = flnm[last(findfirst("n_microbes=",flnm))+1]
    if k == '2'
        k = 22
    elseif k == '8'
        k = 8
    else
        print("bad switch")
    end
    n = flnm[last(findfirst("samplesize=",flnm))+1] == '5' ? flnm[last(findfirst("samplesize=",flnm))+1:last(findfirst("samplesize=",flnm))+3] : flnm[last(findfirst("samplesize=",flnm))+1:last(findfirst("samplesize=",flnm))+4] 

    ix1 = last(findfirst("type=",flnm))+1
    ix2 = first(findfirst(".csv",flnm))-1
    type=flnm[ix1:ix2]
    out=type

    df_tmp = DataFrame(CSV.File("results/simulation/local/realistic-results/$flnm"))

    append!(real_psrf,DataFrame(mu=mu,pi=pi_,R=R,k=k,n=n,type=type,xi=df_tmp[1,:max_xi],gamma=df_tmp[1,:max_gamma]))
end

println("")

display(filter(row -> (row.gamma > 1.2 || row.xi > 1.2), unreal_psrf))

println("")

display(filter(row -> (row.gamma > 1.2 || row.xi > 1.2), real_psrf)) 

println("")