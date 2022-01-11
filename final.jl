using SymPy

# Question (g)
β = 0.995; y_h = 1 ; y_l = 0.75; z_h = 0.3; z_l = 0; ϕ = 0.55; s = 0.02; c = 1; α = 0.45; p = 0.35

θ, w_h, w_l, j_h, j_l = symbols("Θ, w_h, w_l, j_h, j_l")
eq1 = ϕ * y_h + (1-ϕ) * z_h + β* p * ϕ * j_h - w_h 
eq2 = ϕ * y_l + (1-ϕ) * z_l + β* p * ϕ * j_l - w_l 
eq3 = (y_h - w_h) / (1 - β*(1-s)) - j_h
eq4 = (y_l - w_l) / (1 - β*(1-s)) - j_l
eq5 = β * (p / θ) * (0.5 * j_h + 0.5 * j_l) - c
solution = nonlinsolve([eq1, eq2, eq3, eq4, eq5], (θ, w_h, w_l, j_h, j_l))
print(solution)
# θ = 0.524939358937337 and A = p / θ^(1-α)
print("A is ", p / (0.524939358937337^(1-α)))

# Question (h)
print("stationary unemployment is ",s/(0.35+s))

#Question (i)
β = 0.995; y_h = 1 ; y_l = 0.75; z_h = 0.3; z_l = 0; ϕ = 0.55; s = 0.02; c = 1; α = 0.45; A= 0.25

θ, w_h, w_l, j_h, j_l = symbols("Θ, w_h, w_l, j_h, j_l", float=true)
eq1 = ϕ * y_h + (1-ϕ) * z_h + β* A * θ^(1-α) * ϕ * j_h - w_h 
eq2 = ϕ * y_l + (1-ϕ) * z_l + β* A * θ^(1-α) * ϕ * j_l - w_l 
eq3 = (y_h - w_h) / (1 - β*(1-s)) - j_h
eq4 = (y_l - w_l) / (1 - β*(1-s)) - j_l
eq5 = β * (A * θ^(-α)) * (0.5 * j_h + 0.5 * j_l) - c
solution = nonlinsolve([eq1, eq2, eq3, eq4, eq5], (θ, w_h, w_l, j_h, j_l))
print(solution)
# θ = 0.524939358937337 and A = p / θ^(1-α)
print("A is ", p / (0.524939358937337^(1-α)))

