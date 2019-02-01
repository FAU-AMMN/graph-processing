module AbstractVarGraph
#--------------------------------------------------------------------
struct ProblemData
  nType::String
  dimDomain::Integer
  dimRange::Integer
  f
end

struct WeightFunction
  sigma::Float64
  fct::Function
end

struct DistanceFunction
  innerNorm::Function
  outerNorm::Function
  padMethod::String
end

struct Neighborhood
  nType::String
  direction::String
  searchRadius::Integer
  patchRadius::Float64
end
end