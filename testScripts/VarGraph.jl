module VarGraph
using SimpleWeightedGraphs

export constructGraph

include("GridGraph.jl")
using .GridGraph




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function constructGraph(data::Dict, ngh::Dict, distFct::Dict, wFct::Dict)
  if lowercase(ngh["nType"]) == "grid" && lowercase(data["nType"]) == "grid"
    #----------------------------------------------------------------
    #Get the indices that correspond to the neighbors specified by ngh and data in the stencil
    relNghInd = getRelNghInd(lowercase(ngh["direction"]), data["dimDomain"], ngh["searchRadius"])
    #----------------------------------------------------------------
    #Define useful values
    numOfVerts = div(length(data["f"]),data["dimRange"])
    numOfNgh = length(relNghInd)
    numOfEdges = numOfVerts * numOfNgh




    g = SimpleWeightedDiGraph(numOfVerts)
  else
    error("Unsupported neighborhood type!")
  end
  return g
end

end
