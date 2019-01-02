using ProgressMeter
using FFTW, LinearAlgebra
using Plots
gr()
#pyplot()

include("types.jl")
#include("parameters.jl")
#include("times.jl")
#include("mesh.jl")
include("solvers/RK4.jl")
include("initialdata/BellCurve.jl")
include("models/CGBSW.jl")
include("models/Matsuno.jl")
#include("problem.jl")
include("solve.jl")
include("fig.jl")
