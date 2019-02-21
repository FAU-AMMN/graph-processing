using VariationalGraphs
using PyPlot

display("Finished loading!")

dim_domain = 2
num_pts = 1000
const points = rand(dim_domain, num_pts)
const F = rand(num_pts)

include("test_vis.jl")