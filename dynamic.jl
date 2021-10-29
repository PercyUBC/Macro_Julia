using Plots

α = 0.3 
δ = 0.1 
β = 0.9 

kupper = 2 
klower = 0.001 
n = 10 
kgrid = collect(range(klower, stop = kupper, length = n))

v = zeros(n)

tolerance = 0.001
itermax = 20000
iter = 1 
vnew = v .+ 2 * tolerance 

while maximum(abs.(v -vnew)) >  tolerance && iter < itermax
    v = vnew
    c = zeros(n,n)
    for i in 1:n
        for j in 1:n 
            c[i, j] = kgrid[i]^α + (1-δ) * kgrid[i] - kgrid[j]
            if c[i, j] < 0
                c[i, j] = 0
            end 
        end
    end 
    (vnew, cartesianindex) = findmax(log.(c) .+ β * v' ,  dims =2)
    iter = iter +1
end 

kprimeindex = getindex.(cartesianindex, 2)
kprime = kgrid[kprimeindex]

scatter(kgrid, vnew, title="v(k)")






# heatmap(kgrid, kgrid, c, dims = 2, xlab = "k prime", ylab="k",xmirror = true, 
            #yflip=true)