using VariationalGraphs
using Test


dim_domain = 2
num_pts = 10
const points = rand(dim_domain, num_pts)
const F = rand(num_pts)
#------------------------------------------------

include("test_knn.jl")
include("test_ball.jl")