using Plots
using BenchmarkTools

function vfsolve(vnew, kgrid, tolerance, imax)
    α = 0.37
    δ = 0.1
    β = 0.95
    σ = 1.5

    v = vnew .+ 2*tolerance
    cartesianindex = Array{CartesianIndex{2}, length(v)}
    i = 1
    
    c = kgrid.^α + (1-δ)*kgrid .- kgrid'
    c[c .< 0] .= 0
    u_c(x) = (x^(1-σ)- 1)/(1-σ)
    u = u_c.(c)

    while maximum(abs.(v - vnew)) > tolerance && i <= imax
        v = vnew;
        
        (vnew, cartesianindex) = findmax(u .+ β*v', dims = 2);
        i += 1;
    end
    kprimeindex = getindex.(cartesianindex, 2)
    return (v = vnew, kprime = kgrid[kprimeindex], kprimeindex = kprimeindex)
end

# Q1a
kupper = 5
klower = 0.05
n1 = 100
kgrid1 = collect(range(klower, stop = kupper, length = n1))
n2 = 500
kgrid2 = collect(range(klower, stop = kupper, length = n2))
n3 = 100
kgrid3 = collect(range(klower, stop = kupper, length = n3))

(v1, kprime1, kprimeindex1) = vfsolve(zeros(n1), kgrid1, 0.001, 1000);
(v2, kprime2, kprimeindex2) = vfsolve(zeros(n2), kgrid2, 0.001, 1000);
(v3, kprime3, kprimeindex3) = vfsolve(zeros(n3), kgrid3, 0.001, 1000);
plot(kgrid1, v1, label = "v1")
plot!(kgrid2, v2, label = "v2")
plot!(kgrid3, v3, label = "v3")

# Q1b
n=50
kupper = 5
klower = 0.05
kgrid = collect(range(klower, stop=kupper, length=n))
(v1, kprime1, kprimeindex1) = vfsolve(zeros(n), kgrid, 0.0000000001, 1);
(v2, kprime2, kprimeindex2) = vfsolve(zeros(n), kgrid, 0.0000000001, 5);
(v3, kprime3, kprimeindex3) = vfsolve(zeros(n), kgrid, 0.0000000001, 10);

plot(kgrid, v1, label = "v1")
plot!(kgrid, v2, label = "v2")
plot!(kgrid, v3, label = "v3")

function vfsolve(vnew, kgrid, tolerance, imax)
    α = 0.37
    δ = 0.1
    β = 0.95
    σ = 1.5

    v = vnew .+ 2*tolerance
    cartesianindex = Array{CartesianIndex{2}, length(v)}
    i = 1
    
    c = kgrid.^α + (1-δ)*kgrid .- kgrid'
    c[c .< 0] .= 0
    u_c(x) = (x^(1-σ)- 1)/(1-σ)
    u = u_c.(c)

    while maximum(abs.(v - vnew)) > tolerance && i <= imax
        v = vnew;
        
        (vnew, cartesianindex) = findmax(u .+ β*v', dims = 2);
        i += 1;
    end
    kprimeindex = getindex.(cartesianindex, 2)
    print(i)
    return (v = vnew, kprime = kgrid[kprimeindex], kprimeindex = kprimeindex, i=i)
end
vfsolve(zeros(n), kgrid, 0.0000000001, 100000);

function vfsolve(vnew, kgrid, tolerance, imax)
    α = 0.37
    δ = 0.1
    β = 0.95
    σ = 1.5

    v = vnew .+ 2*tolerance
    cartesianindex = Array{CartesianIndex{2}, length(v)}
    i = 1
    
    c = kgrid.^α + (1-δ)*kgrid .- kgrid'
    c[c .< 0] .= 0
    u_c(x) = (x^(1-σ)- 1)/(1-σ)
    u = u_c.(c)

    while maximum(abs.(v - vnew)) > tolerance && i <= imax
        v = vnew;
        
        (vnew, cartesianindex) = findmax(u .+ β*v', dims = 2);
        i += 1;
    end
    kprimeindex = getindex.(cartesianindex, 2)
    return (v = vnew, kprime = kgrid[kprimeindex], kprimeindex = kprimeindex)
end
# 410
# Q1c

n=100
klower = 0.05
kupper = 5
kgrid = collect(range(klower, stop=kupper, length=n))
(v, kprime, kprimeindex) = vfsolve(zeros(n), kgrid, 0.0000000001, 1000);

function k_prime_find(x, grid=kgrid,prime=kprime)
    index = argmin(abs.(x.-grid))
    return prime[index] 
end 
# initialize value 
k_0 = 0.5 
T= 30
Tgrid=collect(1:(T+1))
path=[k_0]
for i in 1:T
    a = k_prime_find(last(path))
    append!(path, a)
end 

scatter(Tgrid, path, label="path of k_t")

#Q1d 
# initialize value 
k_0 = 4.5 
T= 100
Tgrid=collect(1:(T+1))
path=[k_0]
for i in 1:T
    a = k_prime_find(last(path))
    append!(path, a)
end 

scatter(Tgrid, path, label="path of k_t")

#Q1e
function vfsolve(vnew, kgrid, tolerance, imax)
    α = 0.37
    δ = 0.1
    β = 0.95
    σ = 0.5

    v = vnew .+ 2*tolerance
    cartesianindex = Array{CartesianIndex{2}, length(v)}
    i = 1
    
    c = kgrid.^α + (1-δ)*kgrid .- kgrid'
    c[c .< 0] .= 0
    u_c(x) = (x^(1-σ)- 1)/(1-σ)
    u = u_c.(c)

    while maximum(abs.(v - vnew)) > tolerance && i <= imax
        v = vnew;
        
        (vnew, cartesianindex) = findmax(u .+ β*v', dims = 2);
        i += 1;
    end
    kprimeindex = getindex.(cartesianindex, 2)
    return (v = vnew, kprime = kgrid[kprimeindex], kprimeindex = kprimeindex)
end

n=1000
klower = 0.05
kupper = 5
kgrid = collect(range(klower, stop=kupper, length=n))
(v, kprime, kprimeindex) = vfsolve(zeros(n), kgrid, 0.0000000001, 1000);

function k_prime_find(x, grid=kgrid,prime=kprime)
    index = argmin(abs.(x.-grid))
    return prime[index] 
end 

k_0 = 0.5 
T= 30
Tgrid=collect(1:(T+1))
path=[k_0]
for i in 1:T
    a = k_prime_find(last(path))
    append!(path, a)
end 

plot(Tgrid, path, label="path of k_t")

# initialize value 
k_0 = 4.5 
path=[k_0]
for i in 1:T
    a = k_prime_find(last(path))
    append!(path, a)
end 

plot(Tgrid, path, label="path of k_t")
