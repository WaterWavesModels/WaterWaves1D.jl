using ProgressMeter
using FFTW, LinearAlgebra
using Statistics #only for function mean...
ENV["PLOTS_USE_ATOM_PLOTPANE"] = "false"
# Note that javascript-based libraries (for example: PlotlyJS) cannot be shown in the PlotPane due to issues within Atom's internals.
using Plots
gr()
#pyplot()
#plotlyjs()

include("types.jl")
include("solvers/RK4.jl")
include("initialdata/BellCurve.jl")
include("initialdata/HighFreq.jl")
include("initialdata/Random.jl")
include("models/CGBSW.jl")
include("models/CGBSW_naive.jl")
include("models/Matsuno.jl")
include("models/Matsuno_naive.jl")
include("models/Matsuno_mod_naive.jl")
include("models/Boussinesq.jl")
include("models/fdBoussinesq_1.jl")
include("models/fdBoussinesq_1b.jl")
include("models/fdBoussinesq_2.jl")
include("models/fdBoussinesq_2b.jl")
include("models/WaterWaves.jl")
include("LoadSave.jl")
include("fig.jl")
