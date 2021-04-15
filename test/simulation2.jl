using Distributions,Random,TickTock
include("../src/BayesNet.jl")

η  = 1.01
ζ  = 1
ι  = 1
R  = 5
aΔ = 0
bΔ = 0
V = 20
q = floor(Int,V*(V-1)/2)

n=70
π₂=0.3

if π₂==0.8
    println("Case 2")
elseif π₂==0.3
    println("Case 1")
else
    println("Unknown pi value")
end
A  = [zeros(V,V)]
y = zeros(n)

ξ⁰ = map(l -> rand(Bernoulli(π₂)),1:V)
B₀ = zeros(V,V)

for k = 1:V
    for l = (k+1):V
        if ξ⁰[k] == 1 && ξ⁰[l] == 1
            B₀[k,l] = rand(Normal(0.8,1))
            B₀[l,k] = B₀[k,l]
        else
            B₀[k,l] = 0
            B₀[l,k] = 0
        end
    end
end

for i = 1:n
    for k=1:V
        for l = (k+1):V
            A[i][k,l] = rand(Normal(0,1))
            A[i][l,k] = A[i][k,l]
        end
    end
    τ₀² = 1
    ϵᵢ = rand(Normal(0,τ₀²))
    #y[i] = tr(transpose(B₀) * A[i]) + ϵᵢ
    y[i] = sum(upper_triangle(B₀) .* upper_triangle(A[i])) + ϵᵢ
    if i != n
        append!(A,[zeros(V,V)])
    end
end

#nburn = 30000
#nsamp = 20000

nburn = 15000
nsamp = 10000

tick()
τ², u, ξ, γ, D, θ, Δ, M, μ, Λ, πᵥ = BayesNet(A, y, R, nburn=nburn,nsamples=nsamp, V_in = 20, aΔ=1, bΔ=1,ν=10,ι=1,ζ=1)
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

γ_n2 = mean(γ[nburn+1:5:nburn+nsamp])
γ₀ = upper_triangle(B₀)
MSE = 0
for i in 1:190
    global MSE = MSE + (γ_n2[i] - γ₀[i])^2 
end

println("MSE")
println(MSE * (2/(V*(V-1))))
println("")
#show(stdout,"text/plain",γ_n2)


γ_df = DataFrame(n = collect(1:190),l = low, h = high)

sort_df = sort(γ_df,[:l])

#show(stdout,"text/plain",sort_df)
#println("")



show(DataFrame(hat=mean(ξ[nburn:5:nburn+nsamp]),real=ξ⁰))
#show(DataFrame(hat=mean(ξ[nburn:nburn+nsamp]),real=ξ⁰))
println("")
#println(ξ[2])