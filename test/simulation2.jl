using Distributions,Random,TickTock,CSV
include("../src/BayesNet.jl")

η  = 1.01
ζ  = 1
ι  = 1
R  = 5
aΔ = 1
bΔ = 1
ν = 10
V = 20
q = floor(Int,V*(V-1)/2)

data_in = CSV.File("data/simulation2_case2.csv")

X = data_in[:,1:190]
y = data_in[:,191]

#nburn = 30000
#nsamp = 20000

nburn = 15000
nsamp = 10000

tick()
τ², u, ξ, γ, D, θ, Δ, M, μ, Λ, πᵥ = BayesNet(X, y, R, nburn=nburn,nsamples=nsamp, V_in = 20, aΔ=1, bΔ=1,ν=10,ι=1,ζ=1)
tock()


low = zeros(190)
high = zeros(190)
lw = convert(Int64, round(nsamp * 0.1))
hi = convert(Int64, round(nsamp * 0.9))

γ_n = hcat(γ...)

for i in 1:190
    srtd = sort(γ_n[i,nburn+1:nburn+nsamp])
    low[i] = srtd[lw]
    high[i] = srtd[hi]
end

γ_n2 = mean(γ[nburn+1:10:nburn+nsamp])
γ₀ = upper_triangle(B₀)*2
MSE = 0
for i in 1:190
    global MSE = MSE + (γ_n2[i] - γ₀[i])^2
end

println("MSE")
println(MSE * (2/(V*(V-1))))
println("")

println("Other MSE")
println(mean((γ_n2 - γ₀)^2))
#show(stdout,"text/plain",γ_n2)
γ_diff = γ_n2[:,1]-γ₀

show(stdout,"text/plain",DataFrame(gamhat=γ_n2[:,1], gamnaut=γ₀, diff=γ_diff))
#println(typeof(γ_n2))
#println("")
#println(typeof(γ₀))
println("")
#println(typeof(γ_n2-γ₀))

γ_df = DataFrame(n = collect(1:190),l = low, h = high)

sort_df = sort(γ_df,[:l])

#show(stdout,"text/plain",sort_df)
#println("")



show(stdout,"text/plain",DataFrame(hat=mean(ξ[nburn+1:10:nburn+nsamp]),real=ξ⁰))
#show(DataFrame(hat=mean(ξ[nburn:nburn+nsamp]),real=ξ⁰))
println("")
#println(ξ[2])
