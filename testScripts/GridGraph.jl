module GridGraph

using ImageFiltering: Pad, padarray

export getRelNghInd, computePatchDistance2D

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function integralImage(A::AbstractArray)
  intImage = zeros(size(A,1) + 1, size(A,2) + 1)
  intImage[2:end, 2:end] = cumsum(cumsum(A,dims = 1), dims = 2)
  return intImage
end
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
"""
COMPUTEINDEXFROMOFFSET Computes neighbor indices and subscripts from 
 reference offsets.
 This function computes for a given set of offset vectors the
 corresponding subscripts and linear indices with respect to the specified
 padding boundary conditions.

   Input: offsets     - matrix containing the relative neighbor offsets 
                        from their respective reference grid points;
          data        - data struct; should be defined in graph setting script
          k           - number of neighbors per grid point
          pad_method  - the user defined boundary conditions for padding

   Output: subscripts - matrix containing the absolute neighbor subscripts
           indices    - matrix containing the absolute neighbor indices
"""
function computeIndexFromOffset(offsets::AbstractArray, data::Dict, k::Integer, pad_method::String)
  #
  if padMethod != "replicate"
    error("Unsupported padMethod!")
  end
  # extract variables for convenience
  nx = data["nx"];
  ny = data["ny"];
  # check if grid data is three dimensional
  if data["dimDomain"] == 3
    # extract variables for convenience
    nz = data["nz"];
    ndims = data["dimDomain"];
  else # no third dimension
    #set variables for convenience
    nz = 1;
    ndims = 2;
  end
  # generate a meshgrid according to dimensions of data domain
end
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



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
function computePatchDistance2D(data1::Dict, data2::Dict, ngh::Dict,
                                wFct::Dict, dFct::Dict )
  #
  if ngh["padMethod"] != "replicate"
    error("Unsupported padMethod!")
  end
  # Input
  sR = ngh["searchRadius"]
  pR = ngh["patchRadius"]
  innerNorm = dFct["innerNorm"]
  outerNorm = dFct["outerNorm"]
  #
  k = 0;
  max_size = (2*sR + 1)^data1["dimDomain"];
  padding = sR + pR;
  # initialize containers to save distances and offsets of neighbors
  distances = Array{Float64}(undef, data1["ny"], data1["nx"], max_size + k, 1)
  offsets = Array{Float64}(undef, data1["ny"], data1["nx"], max_size + k, 2)
  # 
  f1 = padarray(data1["f"], Pad(:replicate,pR,pR))
  f2 = padarray(data2["f"], Pad(:replicate,padding,padding))
  # left top corner of img
  ltX = ltY = pR + 1 + 1
  #right bottom corner of Image
  rbX = ltX + data1["nx"] - 1
  rbY = ltY + data1["ny"] - 1
  #------------------------------------------------------------------
  # Loop over all shifts in x and y direction
  for x_shift = - sR:sR
    it_x = x_shift + sR + 1
    # left and right index in data2
    data2_lX = -padding + 2 + x_shift
    data2_rX = pR + data2["nx"] + x_shift
    #----------------------------------------------------------------
    for y_shift = -sR:sR
      it_y = y_shift + sR + 1
      # top and bottom index in data2
      data2_tY = -padding + 2 + y_shift
      data2_bY = pR + data2["ny"] + y_shift
      #Restircted data1
      f2_restricted = f2[data2_tY:data2_bY, data2_lX:data2_rX]
      diffsMatrix = innerNorm(f2_restricted - parent(f1))
      integral_images = integralImage(diffsMatrix)
      #
      yRange = ((ltY + pR) : (rbY + pR)) 
      yRange2 = ((ltY + pR - 2 * pR - 1) : (rbY+ pR - 2 * pR - 1)) 
      xRange = ((ltX + pR) : (rbX + pR)) 
      xRange2 = ((ltX + pR - 2 * pR - 1) : (rbX + pR - 2 * pR - 1))
      #
      intImgDiffSums = 
        integral_images[yRange, xRange]# lower right corner
      - integral_images[yRange2, xRange]# upper right corner
      - integral_images[yRange, xRange2]# lower left corner
      + integral_images[yRange2, xRange2]# upper left corner
      #
      intImgDiffSums = outerNorm(intImgDiffSums);
      #
      iteration_mod = mod((it_x - 1) * (sR * 2 + 1) + it_y - 1, max_size) + 1
      #
      for i = 1:data1["ny"], j = 1:data1["nx"]
        offsets[i,j,iteration_mod + k,1] = y_shift
        offsets[i,j,iteration_mod + k,2] = x_shift
      end
      # save distances
      distances[:,:,iteration_mod + k] = intImgDiffSums
    end
  end
  return indices, distances
end
#End computePatchDistance2D
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
