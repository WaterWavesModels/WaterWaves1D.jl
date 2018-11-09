module DeepWaterModels

using Plots
using LinearAlgebra
using FFTW

include("parameters.jl")
include("times.jl")
include("mesh.jl")

export AbstractModel

abstract type AbstractModel end

include("models/cheng.jl")
include("models/matsuno.jl")

include("solvers.jl")
include("fig.jl")

end 
