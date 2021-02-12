using Random,Statistics,Plots
include("../../../Distributions.jl/src/Distributions.jl")

function plot_dists()
    d = Distributions.GeneralizedInverseGaussian(1,2,1.1)
    d2 = Distributions.GeneralizedInverseGaussian(0.9,0.9,0.5)
    d3 = Distributions.GeneralizedInverseGaussian(0.2,0.4,0.4)
    d4 = Distributions.GeneralizedInverseGaussian(2,3,2)

    draws = zeros(100000)
    draws2 = zeros(100000)
    draws3 = zeros(100000)
    draws4 = zeros(100000)
    pdf = zeros(size(0.1:0.25:20))
    pdf2 = zeros(size(0.1:0.25:20))
    pdf3 = zeros(size(0.1:0.25:20))
    pdf4 = zeros(size(0.1:0.25:20))
    for i in 1:100000
        draws[i] = rand(d)
        draws2[i] = rand(d2)
        draws3[i] = rand(d3)
        draws4[i] = rand(d4)
    end
    cnt = 1
    for i in 0.1:0.25:20
        pdf[cnt] = Distributions.pdf(d,i)
        pdf2[cnt] = Distributions.pdf(d2,i)
        pdf3[cnt] = Distributions.pdf(d3,i)
        pdf4[cnt] = Distributions.pdf(d4,i)
        cnt = cnt + 1
    end
    plt = histogram(draws, label="draws", normalize = :pdf, title="a=1, b=2, p=1")
    plot!(plt, 0.1:0.25:20, pdf, label="pdf")
    savefig(plt, "plots/a1b2p1")

    plt2 = histogram(draws2, label="draws", normalize = :pdf, title="a=0.9, b=0.9, p=0.5")
    plot!(plt2, 0.1:0.25:20, pdf2, label="pdf")
    savefig(plt2, "plots/a9b9p5")

    plt3 = histogram(draws3, label="draws", normalize = :pdf, title="a=0.2, b=0.4, p=0.4")
    plot!(plt3, 0.1:0.25:20, pdf3, label="pdf")
    savefig(plt3, "plots/a2b4p4")
    
    plt4 = histogram(draws4, label="draws", normalize = :pdf, title="a=2, b=3, p=2")
    plot!(plt4, 0.1:0.25:20, pdf4, label="pdf")
    savefig(plt4, "plots/a2b3p2")
end

plot_dists()
