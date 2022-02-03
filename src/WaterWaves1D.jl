module WaterWaves1D

# external modules
using ProgressMeter,FFTW
using LinearMaps,IterativeSolvers,LinearAlgebra
using HDF5
import Base.show

# abstract types
export AbstractModel
abstract type AbstractModel end
show(io::IO, m::AbstractModel) =
    try print(io,m.info)
    catch
        print(io,"Model: $(m.label)")
    end

export InitialData
abstract type InitialData end
show(io::IO, i::InitialData) =
    try print(io,i.info)
    catch
        print(io,"Initial data: $(i.label)")
    end

include("solvers/common.jl")

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

include("models/common.jl")
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
