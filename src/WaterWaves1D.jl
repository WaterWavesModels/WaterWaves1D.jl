module WaterWaves1D

# external modules
using ProgressMeter,FFTW
using LinearMaps,IterativeSolvers,LinearAlgebra
using HDF5

# abstract types
export AbstractModel,TimeSolver,InitialData
abstract type AbstractModel end
abstract type TimeSolver end
abstract type InitialData end

# sructures
include("data.jl")
include("mesh.jl")
include("times.jl")
include("init.jl")
include("problem.jl")

# tools
include("tools.jl")
include("figures.jl")

# initial data, models, solvers
include("initialdata/Random.jl")
include("initialdata/CnoidalWaveSerreGreenNaghdi.jl")
include("initialdata/SolitaryWaveSerreGreenNaghdi.jl")
include("initialdata/SolitaryWaveWhithamGreenNaghdi.jl")
include("initialdata/SolitaryWaveWhithamBoussinesq.jl")
include("initialdata/SolitaryWaveWhitham.jl")

include("models/Boussinesq.jl")
include("models/DeepQuadratic.jl")
include("models/IsobeKakinuma.jl")
include("models/Matsuno.jl")
include("models/modifiedMatsuno.jl")
include("models/NonHydrostatic.jl")
include("models/SerreGreenNaghdi.jl")
include("models/SquareRootDepth.jl")
include("models/WaterWaves.jl")
include("models/WhithamBoussinesq.jl")
include("models/WhithamGreenNaghdi.jl")
include("models/WWn.jl")

include("solvers/Euler.jl")
include("solvers/EulerSymp.jl")
include("solvers/RK4.jl")


end
