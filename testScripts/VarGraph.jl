module VarGraph
using SimpleWeightedGraphs

export constructGraph
include("AbstractVarGraph.jl")
using Main.AbstractVarGraph: ProblemData, WeightFunction, DistanceFunction, Neighborhood

include("GridGraph.jl")
using .GridGraph




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function constructGraph(data::ProblemData, ngh::Neighborhood, distFct::DistanceFunction, wFct::WeightFunction)
  if lowercase(ngh.nType) == "grid" && lowercase(data.nType) == "grid"
    #----------------------------------------------------------------
    #Get the indices that correspond to the neighbors specified by ngh and data in the stencil
    relNghInd = getRelNghInd(lowercase(ngh.direction), data.dimDomain, ngh.searchRadius)
    #----------------------------------------------------------------
    #Define useful values
    numOfVerts = div(length(data.f),data.dimRange)
    numOfNgh = length(relNghInd)
    numOfEdges = numOfVerts * numOfNgh




    g = SimpleWeightedDiGraph(numOfVerts)
  else
    error("Unsupported neighborhood type!")
  end
  return g
end

end
