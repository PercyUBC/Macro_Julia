using Plots 
n= 100 
ϵ = rand(n)
plot(1:n, ϵ)
typeof(ϵ)

# for loops
n=100 
ϵ= zeros(n)
for i in 1:n
    ϵ[i]=randn()
end 
plot(1:n, ϵ)

#%%
n = 100 
ϵ = zeros(n)
for i in eachindex(ϵ)
    ϵ[i]= rand()
end 
plot(ϵ)
ϵ[1:5]

ϵ_sum = 0.0 
m = 5 
for ϵ_val in  ϵ[1:m]
    ϵ_sum = ϵ_sum + ϵ_val 
end 
ϵ_mean = ϵ_sum / m 





# functions 
function generating(n)
    ϵ = zeros(n)
    for i in eachindex(ϵ)
        ϵ[i]= (randn())^2 
    end 
    return ϵ
end 

data = generating(10)
plot(data)


f(x)  = x^2 
generation(n)=f.(randn(n))
data = generation(10)

generation(n,fun) = fun.(randn(n))
generation(10, f)

n= 100 
f(x) = x^2 
x = randn(n)
plot(f.(x), label="x^2") 
plot!(x, label = "x")

using Distributions 
function plotthistogram(distribution, n)
    ϵ = rand(distribution, n)
    histogram(ϵ)
end 

lp = Laplace()
plotthistogram(lp, 500)

# while loops 
#find fix point 
using LinearAlgebra 
p = 1 
β = 0.9 
f(x) = p + β * x 
maxiter = 1000
tolerance = 1.0e-8 
initial = 0.3 

nordiff = Inf 
oldv= initial 
iter = 1 
while nordiff > tolerance && iter <= maxiter
    newv= f(oldv)
    nordiff = norm(newv - oldv)
    oldv = newv 
    iter = iter + 1 
end 
println("Fixed point is $oldv, and |f(x)-x|=$nordiff.")

# for style 
p = 1 
β = 0.9 
f(x) = p + β * x 
maxiter = 1000
tolerance = 1.0e-8 
initial = 0.3 

nordiff = Inf 
oldv= initial 
iter = 1 
for i in 1:maxiter
    newv=f(oldv)
    nordiff = norm(newv - oldv)
    if nordiff < tolerance
        iter = i 
        break 
    end 
    oldv = newv 
end
println("Fixed point is $oldv, and |f(x)-x|=$nordiff in $iter iterations")

# better way 
function fixpoint(maxiter, tolerance, initial, p, β)
    f(x) = p + β * x
    iter = 1 
    v_old = initial
    normdiff = Inf
    while normdiff > tolerance && iter <= maxiter
        v_new = f(v_old)
        normdiff = norm(v_new - v_old)
        v_old = v_new 
        iter += 1 
    end 
    return (v_old, normdiff, iter)
end 
v_old, normdiff, iter = fixpoint(10000, 1.0E-8, 1, 1, 0.9)
println("Fixed point is $v_old, and |f(x)-x|=$normdiff in $iter iterations")

# better style 
p = 1 
β = 0.9
f(x) = p + β * x
function fixpoint_new(maxiter, tolerance, initial, g)
    iter = 1 
    v_old = initial
    normdiff = Inf
    while normdiff > tolerance && iter <= maxiter
        v_new = g(v_old)
        normdiff = norm(v_new - v_old)
        v_old = v_new 
        iter += 1 
    end 
    return (v_old, normdiff, iter)
end
v_old, normdiff, iter = fixpoint_new(10000, 1.0E-8, 1, f)
println("Fixed point is $v_old, and |f(x)-x|=$normdiff in $iter iterations")

# best style 
function fixpoint_final(f; initial, maxiter=10000, tolerance=1.0E-8)
    v_old = initial
    iter = 1 
    normdiff = Inf 
    while normdiff > tolerance && iter <= maxiter
        v_new = f(v_old)
        normdiff = norm(v_new - v_old)
        v_old = v_new 
        iter += 1 
    end 
    return (value=v_old, diff=normdiff, iteration=iter)
end 
r = 2 
f(x) =  r * x * (1 - x)
sol = fixpoint_final(f, initial=0.8)
println("Fixed point is $(sol.value), and |f(x)-x|=$(sol.diff) in $(sol.iteration) iterations")

using NLsolve 
p = 1.0 
β = 0.9 
sol = fixedpoint(v -> p .+ β * v, [1.0])

println("Fixed point is $(sol.zero), and |f(x)-x|=$(norm(sol.zero - f(sol.zero))) in $(sol.iterations) iterations")

p = [1.0, 2.0, 0.1]
β = 0.9 
ini = [1.0, 1.0, 1.0]
sol = fixedpoint(v -> p .+ β * v, ini)
println("Fixed point is $(sol.zero), and |f(x)-x|=$(norm(sol.zero - f(sol.zero))) in $(sol.iterations) iterations")
