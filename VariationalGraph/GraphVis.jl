module GraphVis
using Plots
using GraphRecipes

export vis

function vis(edges::AbstractArray, points::AbstractArray)
    #plt = plot(points[1,:], points[2,:], seriestype=:scatter, legend = false)
    if size(points,1) == 2
        graphplot(edges[:,1], edges[:,2], x = points[1,:], y = points[2,:])
    elseif size(points,1) == 3
        graphplot(edges[:,1], edges[:,2], x = points[1,:], y = points[2,:], z = points[3,:], dim = 3)
    else
        display("Unsupported dimesion!")
    end
end


end