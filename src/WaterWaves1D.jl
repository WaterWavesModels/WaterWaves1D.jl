module WaterWaves1D

# external modules
using ProgressMeter,FFTW
using LinearMaps,IterativeSolvers,LinearAlgebra
using HDF5
import Base.show

# abstract types
export AbstractModel
"""
Abstract type whose subtypes are the models from which initial-value problems can be built,
through `Problem( model :: AbstractModel, initial :: InitialData, param :: NamedTuple )`
"""
abstract type AbstractModel end
show(io::IO, m::AbstractModel) =
    try print(io,m.info)
    catch
        print(io,"Model: $(m.label)")
    end

export InitialData

"""
Abstract type defining initial data from which initial-value problems can be built,
through `Problem( model :: AbstractModel, initial :: InitialData, param :: NamedTuple )`
"""
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
include("recipes.jl")

# initial data, models, solvers
include("initialdata/Random.jl")
include("initialdata/CnoidalWaveSerreGreenNaghdi.jl")
include("initialdata/SolitaryWaveSerreGreenNaghdi.jl")
include("initialdata/SolitaryWaveWhithamGreenNaghdi.jl")
include("initialdata/SolitaryWaveWhithamBoussinesq.jl")
include("initialdata/SolitaryWaveWhitham.jl")

include("models/common.jl")
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

include("solvers/Euler.jl")
include("solvers/EulerSymp.jl")
include("solvers/RK4.jl")


end
