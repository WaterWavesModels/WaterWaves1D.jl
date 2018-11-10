module DeepWaterModels

using Plots
using LinearAlgebra
using FFTW

export AbstractModel
abstract type AbstractModel end

include("times.jl")
include("mesh.jl")
include("solvers.jl")
include("parameters.jl")
include("initial_data.jl")
include("models/cheng.jl")
include("models/matsuno.jl")
include("problem.jl")
include("fig.jl")

end 
