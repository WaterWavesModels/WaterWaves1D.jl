export AbstractModel
export TimeSolver
export InitialData

abstract type AbstractModel end
abstract type TimeSolver end
abstract type InitialData end

include("data.jl")
include("times.jl")
include("mesh.jl")
include("problem.jl")

