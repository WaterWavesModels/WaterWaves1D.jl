module DeepWaterModels

using Plots
using LinearAlgebra
using FFTW
using ProgressMeter

export AbstractModel
abstract type AbstractModel end

include("parameters.jl")
include("times.jl")
include("mesh.jl")
include("solvers.jl")
include("initial_data.jl")
include("models/cheng.jl")
include("models/matsuno.jl")
include("problem.jl")
include("fig.jl")

end 
