# Question 1 
function factorial2(n)
    value = 1
    for i in 1:n
        value = value * i 
    end 
    return value 
end 

println(factorial2(10))
factorial(10) == factorial2(10) && println("The function is perfect.")

# Question 2 
function binomial_rv(n, p)
    value = 0 
    for i in 1:n
        rv = rand()
        if rv<p
            value += 1 
        end 
    end 
    return value 
end 


# Question 3 
iteration = 100000000
frequency=0
for i in 1:iteration
    unit = randn(2)
    norm_length = norm(unit - [0.5, 0.5])
    if norm_length <= 0.5 
        frequency += 1 
    end 
end 
π_estimate = frequency/iteration/ 0.5^2 
println("The π is estimated to be $π_estimate in $iteration iterations")

# Question 4 
function award()
    sum = 0
    realization = rand(10)
    for i in eachindex(realization)
        realization[i] > 0.5 ? sum += 1 : sum = 0 
        if sum >2
            println("You get one dollar!!!") 
            break
        end 
    end
    if sum < 3 
        println("You can try one more time.")
    end 
end 

award()
     
        

