using DelimitedFiles
using VariationalGraphs
using LightGraphs
#------------------------------------------------
display("Finished loading!")
#------------------------------------------------
d = 3
dim_domain = 2
num_pts = N = 50
# const POINTS = rand(dim_domain, num_pts)
const POINTS = readdlm("../data/points.txt", Float64)
#------------------------------------------------
num_nghs = 3
epsilon = 1
#------------------------------------------------
inner_norm(x) = x.^2
outer_norm(x) = sqrt(x)
w_fct(x) = 1 ./ (x.^2)
#------------------------------------------------
ngh = Dict([("nType", "knn"), ("direction", "full"), ("searchRadius", 1), ("num_nghs", num_nghs), ("epsilon", epsilon)])
distFct = Dict([("inner_norm", inner_norm), ("outer_norm", outer_norm)])
wFct = Dict([("sigma", 0), ("fct", w_fct)])
data = Dict([("type", "point_cloud"), ("dimDomain", dim_domain), ("dimRange", 1),("points", POINTS)])
#------------------------------------------------
g = constructGraph(data,ngh,distFct,wFct)
#------------------------------------------------
alpha = 0.001
data = randn(g.num_verts,d)
data = data .- minimum(data)
data = data./maximum(data)
const W = ones(g.num_edges)
#const F = data[:]
tmp = readdlm("../data/F.txt", Float64)
const f0 = reshape(tmp,150)
F = similar(f0)
for i = 1:d
    mean = sum(f0[((i - 1) * N + 1):(i* N)])/N
    for j = ((i - 1) * N + 1):(i* N)
        F[j] = mean
    end
end
#------------------------------------------------
include("test_generate_flowgraph.jl")