#cd("bayesian_network_regression_imp")
using DrWatson;@quickactivate
using Distributed#,ProfileView,Profile
addprocs(2;topology="master_worker")
include("simulation/run_simulation.jl")
simnum = 2
πₛ = 0.8
k = 22
sampsize = 500
μₛ = 0.8
loadinfo = Dict("simnum"=>simnum,"pi"=>πₛ,"mu"=>μₛ,"n_microbes"=>k,"out"=>"xis","samplesize"=>sampsize)
loadinfo["out"] = "XYs"
loadinfo["type"] = "additive_random"
loadinfo["edge_mu"] = 0.4
simtypes = Dict(1 => "unrealistic", 2 => "realistic")
data_in = DataFrame(CSV.File(datadir(joinpath("simulation",simtypes[simnum]),savename(loadinfo,"csv",digits=1))))
loadinfo["out"] = "bs"
b_in = DataFrame(CSV.File(datadir(joinpath("simulation",simtypes[simnum]),savename(loadinfo,"csv",digits=1))))
B₀ = convert(Array{Float64,1},b_in[!,:B])
X = Matrix(data_in[:,names(data_in,Not("y"))])
y = SVector{size(X,1)}(data_in[:,:y])
q = size(X,2)
V = convert(Int,(1 + sqrt(1 + 8*q))/2)
aΔ=1.0;bΔ=1.0;ν=10;R=9;ι=1.0;ζ=1.0;η=1.01
nburn=30000;nsamp=20000

#### in parallel
states = GenerateSamples!(X, y, R, η=η, nburn=nburn,nsamples=nsamp, V = V, aΔ=aΔ, bΔ=bΔ,ν=ν,ι=ι,ζ=ζ,x_transform=false)

#@profview GenerateSamples!(X, y, R, η=η, nburn=nburn,nsamples=nsamp, V = V, aΔ=aΔ, bΔ=bΔ,ν=ν,ι=ι,ζ=ζ,x_transform=false)
#@profile GenerateSamples!(X, y, R, η=η, nburn=nburn,nsamples=nsamp, V = V, aΔ=aΔ, bΔ=bΔ,ν=ν,ι=ι,ζ=ζ,x_transform=false)
#Profile.print()
### in sequence
#GenerateSamples!(X, y, R, η=η, nburn=nburn,nsamples=nsamp, V = V, aΔ=aΔ, bΔ=bΔ,ν=ν,ι=ι,ζ=ζ,x_transform=false,in_seq=true)
#exit()
DataFrame(states.psrf)
exit()
