module Con2Prim

using LinearAlgebra
using NonlinearSolve
using StaticArrays

include("boxes.jl")
include("solve1d.jl")
include("solve2d.jl")
include("solve3d.jl")
include("hydro.jl")

end
