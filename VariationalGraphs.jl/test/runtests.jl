using VariationalGraphs
using DelimitedFiles
using Test


dim_domain = 2
num_pts = 10
const points = rand(dim_domain, num_pts)
const F = rand(num_pts)
tmp = readdlm("2D_PointData", ',', Float64)
const POINTS_2D = collect(tmp')
#------------------------------------------------

include("test_knn.jl")
include("test_ball.jl")