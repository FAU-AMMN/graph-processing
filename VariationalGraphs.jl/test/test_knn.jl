
@testset "basic" begin
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
    data = Dict([("type", "point_cloud"), ("dimDomain", dim_domain), ("dimRange", 1), ("f", F),("points", points)])
    #------------------------------------------------
    G = constructGraph(data,ngh,distFct,wFct)
    #------------------------------------------------ 
    H = undirect(G)
    A = getadjMat(G)
    B = getadjMat(H)
    
    @test B == B'
    @test 0.5*(A + A') == B
end