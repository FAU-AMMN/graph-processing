# VariationalGraphs.jl
Implemantation of some algorithms concerning variational problems on graphs.

## Graph generators
ToDo
## Graph algorithms
### Cut pursuit
Procedure for a cut pursuit algorithm with a primal dual approach
```jl
    cutpursuit(d, f0, g, w, cutAlpha, nIter, cut)
```
* `d:`
A parameter spcifiying the dimension of the given data. This is due to the 
development process of this project, as the data f0 has to be passed as a 
one dimenionsl array, and hence the information is lost during the reshape 
process. TR: We probably drop this procedure.
- typeof(d) = Int64

* `f0:`
The Original input data. At this stage 
multidimensional input has to be rehaped into this format, which might change 
later. typeof(f0) = Array{T, 1} where T <: Real.

* `g:`
This variable stores the underlying graph, the algorithm should be 
performed on. typeof(g) = VariationalGraph{T} where T <: Real.

* `w:`
Weights corresponding to the edges in the variationalGraph g. This is 
currently a bit unclear, as g itself stores weights. It should be 
discussed why g needs to store weights and how the correspond to w. 
typeof(w) = Array{T, 1} where T <: Real.

* `cutAlpha:`
Regularisation parameter to generate the flow graph. 
typeof(cutAlpha) = T where T <: Real. 

* `nIter:`
The Number of outer iterations for the cut pursuit algorithm. 
typeof(nIter) = Int64.

* `cut:`
This parameter allows to dispatch to different choices of cut types for 
the algorithm. Below the possible options are listed. 
typeof(cut) <: cut_type.

* `aniso():`
Perform an anisotropic cut for the cut pusruit algorithm.

* `iso():`
Currently unsupported.

### Primal Dual
