#------------------------------------------------
num_nghs = 3
epsilon = 60 * 1/num_pts
#------------------------------------------------
inner_norm(x) = x.^2
outer_norm(x) = sqrt(x)
w_fct(x) = 1 ./ (x.^2)
#------------------------------------------------
ngh = Dict([("nType", "epsilon"), ("direction", "full"), ("searchRadius", 1), ("num_nghs", num_nghs), ("epsilon", epsilon)])
distFct = Dict([("inner_norm", inner_norm), ("outer_norm", outer_norm)])
wFct = Dict([("sigma", 0), ("fct", w_fct)])
data = Dict([("type", "point_cloud"), ("dimDomain", dim_domain), ("dimRange", 1), ("f", F),("points", points)])
#------------------------------------------------
G = constructGraph(data,ngh,distFct,wFct)
#------------------------------------------------ 
colors = ["b" "r" "c" "m" "k" "g"]
# PyPlot.scatter(points[1, :], points[2, :])

edges = G.edges


@time for u = 1:G.num_verts
    c = colors[rand(1:length(colors))]
    lines = points[:,edges[u]]
    plot(points[1,u], points[2,u], color = c, marker = "o")
    #annotate(string(u), (points[1,u], points[2,u]))
    for i = 1:length(edges[u])
        v = edges[u][i]
        plot((points[1, u], points[1, v]), (points[2, u], points[2, v]), color = c)
    end
end
