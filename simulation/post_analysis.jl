

function calc_msey()
    loadinfo = Dict("simnum"=>1,"pi"=>0.3,"mu"=>0.8,"n_microbes"=>8,"out"=>"xis","samplesize"=>100)
    loadinfo["out"] = "XYs"
    
    data_in = DataFrame(CSV.File(datadir(joinpath("simulation","unrealistic"),savename(loadinfo,"csv",digits=1))))

    loadinfo["out"] = "As"
    A = DataFrame(CSV.File(datadir(joinpath("simulation","unrealistic"),savename(loadinfo,"csv",digits=1))))

    loadinfo["R"] = 9
    loadinfo["nu"] = 10

    loadinfo["out"] = "edges"
    b = DataFrame(CSV.File(projectdir("results","simulation","unrealistic",savename(loadinfo,"csv",digits=1))))

    loadinfo["out"] = "mu"
    μ = DataFrame(CSV.File(projectdir("results","simulation","unrealistic",savename(loadinfo,"csv",digits=1))))

    X = Matrix(data_in[:,names(data_in,Not("y"))])
    y = SVector{size(X,1)}(data_in[:,:y])

    n = size(y,1)
    y_pred = zeros(n)
    μ_post = μ[1,"mean"]
    A_mat = Matrix(A)
    b_a = Vector(b["mean"])
    for i in 1:n
        y_pred[i] = μ_post + sum(b_a .* A_mat[i,:])
    end
    MSEy = sum((y - y_pred).^2)/n
    println(MSEy)
end