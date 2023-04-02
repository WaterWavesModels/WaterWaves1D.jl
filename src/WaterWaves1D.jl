module WaterWaves1D

# external modules
using ProgressMeter,FFTW
using LinearMaps,IterativeSolvers,LinearAlgebra
using HDF5
import Base.show

# abstract types and structures
include("mesh.jl")
include("times.jl")
include("data.jl")
include("init.jl")
include("solver.jl")
include("models.jl")
include("problem.jl")

# tools
include("loadsave.jl")
include("tools.jl")
include("recipes.jl")

# initial data, models, solvers
include("initialdata/Random.jl")
include("initialdata/CnoidalWaveSerreGreenNaghdi.jl")
include("initialdata/SolitaryWaveSerreGreenNaghdi.jl")
include("initialdata/SolitaryWaveWhithamGreenNaghdi.jl")
include("initialdata/SolitaryWaveWhithamBoussinesq.jl")
include("initialdata/SolitaryWaveWhitham.jl")

include("models/Airy.jl")
include("models/Boussinesq.jl")
include("models/AkersNicholls.jl")
include("models/IsobeKakinuma.jl")
include("models/Matsuno.jl")
include("models/modifiedMatsuno.jl")
include("models/NonHydrostatic.jl")
include("models/SerreGreenNaghdi.jl")
include("models/SaintVenant.jl")
include("models/SquareRootDepth.jl")
include("models/WaterWaves.jl")
include("models/WhithamBoussinesq.jl")
include("models/WhithamGreenNaghdi.jl")
include("models/WWn.jl")

include("models/relaxedGreenNaghdi.jl")
include("models/SVN.jl")
include("models/Toy.jl")

include("solvers/Euler.jl")
include("solvers/EulerSymp.jl")
include("solvers/RK4.jl")

end
