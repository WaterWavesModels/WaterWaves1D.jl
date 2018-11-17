module DeepWaterModels

using Plots
using LinearAlgebra
using FFTW
using ProgressMeter

export AbstractModel
export TimeSolver
export InitialData

abstract type AbstractModel end
abstract type TimeSolver end
abstract type InitialData end


include("parameters.jl")
include("times.jl")
include("mesh.jl")
include("solvers/RK4.jl")
include("initialdata/Bump.jl")
include("models/CGBSW.jl")
include("models/Matsuno.jl")
include("problem.jl")
include("fig.jl")

end
