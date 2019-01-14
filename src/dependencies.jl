using ProgressMeter
using FFTW, LinearAlgebra
using Plots
gr()
#pyplot()

include("types.jl")
include("solvers/RK4.jl")
include("initialdata/BellCurve.jl")
include("initialdata/HighFreq.jl")
include("models/CGBSW.jl")
include("models/CGBSW_naive.jl")
include("models/Matsuno.jl")
include("models/Matsuno_naive.jl")
include("models/Matsuno_mod_naive.jl")
include("LoadSave.jl")
include("solve.jl")
include("fig.jl")
