module GridGraph

export getRelNghInd

include("AbstractVarGraph")
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Get assciated neighbors
function getRelNghInd(direction::String, dimDomain::Integer, sR::Integer)
  #------------------------------------------------------------------
  #Direct Neighbors
  if (direction == "directneighbors") && (sR == 1)
    if dimDomain == 2
      relNghInd = [2,4,6,8]
    elseif dimDomain == 3
      relNghInd = [5,11,13,15,17,23]
    else
      error("Unsupported Domain Dimension, must be in {2,3}!")
    end
  #------------------------------------------------------------------
  #Forward Neighbors
  elseif (direction == "forward") && (sR == 1)
    if dimDomain == 2
      relNghInd = [6,8]
    elseif dimDomain == 3
      relNghInd = [15,17,23]
    else
      error("Unsupported Domain Dimension, must be in {2,3}!")
    end
  #------------------------------------------------------------------
  #Backward Neighbors
  elseif (direction == "backward") && (sR == 1)
    if dimDomain == 2
      relNghInd = [2,4]
    elseif dimDomain == 3
      relNghInd = [5,11,13]
    else
      error("Unsupported Domain Dimension, must be in {2,3}!")
    end
  #------------------------------------------------------------------
  #Full
  elseif (direction == "full")
    indices = collect(1:(2*sR + 1)^dimDomain)
    center = div((2*sR + 1)^dimDomain,2)
    relNghInd = vcat(collect(1:(center - 1)),collect((center + 1) : (2*sR + 1)^dimDomain))
  #------------------------------------------------------------------
  #Unknown
  else
    error("Unsupported neighborhood configuration!")
  end
  return relNghInd
end
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Code
function computePatchDistance2D(data1::ProblemData, data2::ProblemData, ngh::Neighborhood, 
                                wFct::WeightFunction, dFct::DistanceFunction )
  sR = ngh.searchRadius
  pR = ngh.patchRadius
  
  k = 0;
  max_size = (2*sR + 1)^data1.dimDomain;
  padding = sR + pR;
end
end
