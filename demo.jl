my_integer = 5
my_float = 5.
my_range = 1:7
my_vector = [1, 2, 3]
my_string = "hello"
my_character = 'o'
my_boolean = true
println("Julia is good")
println("It's even better than python!")
α = 5 
print("hellow world")
print("hellow world")
print("hellow world")

@show α + my_integer

A = [1, 2, 3]   # column vector (array)
B = [1 2 3]     # row vector (array)
ones(3)
zeros(5)
C = [1 3; 5 9]  # array

A = 1:9
B = 1:0.5:9

typeof(A)
eltype(A)
size(A)
length(A)
ndims(A)

A[1]    # first element
A[end]  # last element
C[1,2]  # element in row 1, column 2
C[3]    # Julia iterates column-wise
C[3,]   # Same as above
C[2,:]  # Select row


using LinearAlgebra
I
C*I

C'
C^2
C.^2
inv(C)
I/C

using Random
using Random
rand(1)
rand(1, 1)
randn(2)
Random.seed!(1234)
rand(3)
Random.seed!(1234)
rand(1)
