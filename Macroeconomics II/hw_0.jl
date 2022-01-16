using Random, Distributions, Statistics, Plots, StatsBase
using StatsBase
using Printf

function tauchen(ρ, σ, N, m=3)
    dist = Normal(0,1)
    kupper = m * σ / (1 - ρ^2)^(1/2)
    klower = - m * σ / (1 - ρ^2)^(1/2)
    state = collect(range(klower, stop=kupper, length=N))
    w = state[2] - state[1]
    function trans_prob(j, k)
        kright = (state[k] - ρ* state[j] + w/2) / σ
        kleft = (state[k] - ρ* state[j] - w/2) / σ
        if (k>1 & k<N) 
            return (cdf(dist, kright) - cdf(dist, kleft))
        elseif k==1 
            return cdf(dist, kright)
        else
            return cdf(dist, kleft)
        end
    end
    trans_matrix = zeros(N, N)
    for a in 1:N
        for b in 1: N
            trans_matrix[a, b]= trans_prob(a,b)
        end 
    end 
    return (state, trans_matrix)
end


function rouwenhorst(ρ, σ, N, m=3)
    kupper = m * σ / (1 - ρ^2)^(1/2)
    klower = - m * σ / (1 - ρ^2)^(1/2)
    state = collect(range(klower, stop=kupper, length=N))
    K = 2
    #p = tauchen(ρ, σ, 2)[2][1,1]
    #q = tauchen(ρ, σ, 2)[2][2,2]
    p = (1+ρ)/2
    q = (1+ρ)/2
    
    trans_matrix = [p (1-p); (1-q) q]
    while K < N
        mtx = zeros(K+1, K+1)
        mtx_1 = copy(mtx)
        mtx_1[1:K, 1:K] = copy(trans_matrix)
        mtx_2 = copy(mtx)
        mtx_2[1:K, 2:K+1] = copy(trans_matrix)
        mtx_3 = copy(mtx)
        mtx_3[2:K+1, 1:K] = copy(trans_matrix)
        mtx_4 = copy(mtx)
        mtx_4[2:K+1, 2:K+1] = copy(trans_matrix)
        trans_matrix = mtx_1 .* p + mtx_2 .* (1-p) + mtx_3 .* (1-q) + mtx_4 .* q
        trans_matrix[2:end-1, :] = trans_matrix[2:end-1, :] .* 0.5
        K = K+1
    end    
    return state, trans_matrix
end


tauchen(0.2, 0.4, 6)[2][2,:]
rouwenhorst(0.5, 1, 6)[2]


function gen_tauchen(ρ, σ, N, t, m=3)
    state, mat = tauchen(ρ, σ, N, m)
    k_0 = (N - N%2)/2 + 1
    y_0 = state[k_0]
    k_path = [k_0]
    y_path=[y_0]
    for i in 2:t
        k_new = wsample(collect(1:N), mat[last(k_path), :], 1)
        append!(k_path, k_new)
        append!(y_path, state[k_new])
    end   
    return y_path  
end


function gen_rouwenhorst(ρ, σ, N, t, m=3)
    state, mat = rouwenhorst(ρ, σ, N, m)
    k_0 = convert(Int64, (N - N%2)/2 + 1)
    y_0 = state[k_0]
    k_path = [k_0]
    y_path=[y_0]
    for i in 2:t
        k_new = wsample(collect(1:N), mat[last(k_path), :], 1)
        append!(k_path, k_new)
        append!(y_path, state[k_new])
    end   
    return y_path  
end

function original(ρ, σ,t)
    y_0 = 2.0
    y_path = zeros(t)
    y_path[1]=y_0
    for i in 2:t
        y_path[i] = y_path[i-1]* ρ + rand(Normal(0, σ))
    end   
    return y_path      
end


function percy_summary(d)
    m = mean(d)
    v = var(d)
    auto=autocor(d, [1], demean=true)[1]
    qant_25 = quantile(d, 0.25)
    qant_75 = quantile(d, 0.75)
    println("Mean is $m")
    println("Variance is $v")
    println("25% quantile is $qant_25")
    println("75% quantile is $qant_75")
    println("Auto-correlation is $auto")
end

# Q3 original simulation
y=original(0.2, 0.4, 1000)
println("For orginial equation simulation")
percy_summary(y)


# Q3 tauchen simulation
y=gen_tauchen(0.2, 0.4, 3, 1000)
println("For Tauchen simulation")
percy_summary(y)

# Q3 rouwenhorst simulation
y=gen_rouwenhorst(0.2, 0.4, 3, 1000)
println("For Rouwenhorst simulation")
percy_summary(y)



# Q4 tauchen simulation
y=gen_tauchen(0.2, 0.4, 3, 1000)
println("For Tauchen simulation")
percy_summary(y)

# Q4 rouwenhorst simulation
y=gen_rouwenhorst(0.2, 0.4, 3, 1000)
println("For Rouwenhorst simulation")
percy_summary(y)


# Q5 
#N=3, ρ=0.7
y=gen_tauchen(0.7, 0.4, 3, 1000)
println("For Tauchen simulation")
percy_summary(y)

y=gen_rouwenhorst(0.7, 0.4, 3, 1000)
println("For Rouwenhorst simulation")
percy_summary(y)

#N=10, ρ=0.7
y=gen_tauchen(0.7, 0.4, 10, 1000)
println("For Tauchen simulation")
percy_summary(y)

y=gen_rouwenhorst(0.7, 0.4, 10, 1000)
println("For Rouwenhorst simulation")
percy_summary(y)

########################################
#N=3, ρ=0.9
y=gen_tauchen(0.9, 0.4, 3, 1000)
println("For Tauchen simulation")
percy_summary(y)

y=gen_rouwenhorst(0.9, 0.4, 3, 1000)
println("For Rouwenhorst simulation")
percy_summary(y)

#N=10, ρ=0.9
y=gen_tauchen(0.9, 0.4, 10, 1000)
println("For Tauchen simulation")
percy_summary(y)

y=gen_rouwenhorst(0.9, 0.4, 10, 1000)
println("For Rouwenhorst simulation")
percy_summary(y)

########################################
#N=3, ρ=0.98
y=gen_tauchen(0.98, 0.4, 3, 1000)
println("For Tauchen simulation")
percy_summary(y)

y=gen_rouwenhorst(0.98, 0.4, 3, 1000)
println("For Rouwenhorst simulation")
percy_summary(y)

#N=10, ρ=0.98
y=gen_tauchen(0.98, 0.4, 10, 1000)
println("For Tauchen simulation")
percy_summary(y)

y=gen_rouwenhorst(0.98, 0.4, 10, 1000)
println("For Rouwenhorst simulation")
percy_summary(y)