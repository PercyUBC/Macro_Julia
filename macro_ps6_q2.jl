using Plots
using BenchmarkTools
function crra(x, σ)
    if σ ==1 
        a = log(x) 
    else 
        a = (x^(1-σ) - 1)/(1-σ) 
    end 
    return a
end 


function vfsolve(vnew, kgrid, tolerance, imax, σ)
    β = 0.9
    v = vnew .+ 2*tolerance
    cartesianindex = Array{CartesianIndex{2}, length(v)}
    i = 1
    
    c = kgrid.-kgrid'
    c[c .< 0 ] .= 0
    u = crra.(c, σ)

    while maximum(abs.(v - vnew)) > tolerance && i <= imax
        v = vnew;
        
        (vnew, cartesianindex) = findmax(u .+ β*v', dims = 2);
        i += 1;
    end
    kprimeindex = getindex.(cartesianindex, 2)
    return (v = vnew, kprime = kgrid[kprimeindex], kprimeindex = kprimeindex)
end

# Q2b
n = 1000 
klower = 0.01
kupper = 5
kgrid = collect(range(klower, stop=kupper, length=n))
(v, kprime, kprimeindex) = vfsolve(zeros(n), kgrid, 0.001, 100000, 1.5)
plot(kgrid, v, label="v(k)")

# Q2c
(v, kprime, kprimeindex) = vfsolve(zeros(n), kgrid, 0.001, 100000, 1.25)
plot(kgrid, v, label="σ=1.25")
(v, kprime, kprimeindex) = vfsolve(zeros(n), kgrid, 0.001, 100000, 1.1)
plot!(kgrid, v, label="σ=1.1")
(v, kprime, kprimeindex) = vfsolve(zeros(n), kgrid, 0.001, 100000, 1.001)
plot!(kgrid, v, label="σ=1.001")

# Q2d
(v, kprime, kprimeindex) = vfsolve(zeros(n), kgrid, 0.001, 100000, 1.5)
function k_prime_find(x, grid=kgrid,prime=kprime)
    index = argmin(abs.(x.-grid))
    return prime[index] 
end 
# initialize value 
k_0 = 5.0 
T= 15
Tgrid=collect(1:(T+1))
path=[k_0]
for i in 1:T
    a = k_prime_find(last(path))
    append!(path, a)
end 

plot(Tgrid, path, label="path of k_t")


