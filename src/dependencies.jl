using ProgressMeter
using FFTW, LinearAlgebra
using LinearMaps,IterativeSolvers
using Statistics #only for function mean in WaterWaves.jl...
using JLD #only for LoadSave
ENV["PLOTS_USE_ATOM_PLOTPANE"] = "false"
# Note that javascript-based libraries (for example: PlotlyJS) cannot be shown in the PlotPane due to issues within Atom's internals.
using Plots
gr()
#pyplot()
#plotlyjs()

include("types.jl")

include("solvers/RK4.jl")
include("solvers/RK4_naive.jl")

include("initialdata/Init.jl")
include("initialdata/BellCurve.jl")
include("initialdata/HighFreq.jl")
include("initialdata/Random.jl")
include("initialdata/SolitaryWaveWhitham.jl")
include("initialdata/SolitaryWaveWhithamBoussinesq.jl")
include("initialdata/SolitaryWaveWhithamGreenNaghdi.jl")
include("initialdata/CnoidalWaveWhithamGreenNaghdi.jl")

include("models/CGBSW.jl")
include("models/CGBSW_naive.jl")
include("models/Matsuno.jl")
include("models/Matsuno_naive.jl")
include("models/Matsuno_mod_naive.jl")
include("models/Boussinesq.jl")
include("models/PseudoSpectral.jl")
include("models/WaterWaves.jl")
include("models/WhithamBoussinesq.jl")
include("models/WhithamGreenNaghdi.jl")
include("models/WhithamGreenNaghdiGPU.jl") #comment this line if you have problems with CuArrays

include("LoadSave.jl")
include("fig.jl")
include("tools.jl")
