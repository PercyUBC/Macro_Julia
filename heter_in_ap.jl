using Plots
using BenchmarkTools

# quadratic function plot
n = 1000 
klower = -2.7
kupper = 3.3
kgrid = collect(range(klower, stop=kupper, length=n))
function quadratic(x)
    index = 1 .+ ( x .+ 0.5) .* (x .- 1.1)
    return index 
end 

curve = quadratic(kgrid)

plot(kgrid, curve, label = "r", xlabel="η", dpi=300)
plot!([0], seriestype="vline", label="", dpi=300)
plot!([0.3], seriestype="vline", label="", dpi=300)