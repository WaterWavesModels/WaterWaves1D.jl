export ProblemSave, save, load
using JLD
struct ProblemSave

    model   :: Symbol
    initial :: Symbol
    param   :: Dict
    solver  :: Symbol
    data    :: Data

    function ProblemSave(model, initial, param, solver, data)

         new( model, initial, param, solver, deepcopy(data))

    end

end

import Base.convert

function convert(::Type{ProblemSave}, p :: Problem )

    model   = Symbol(typeof(p.model))
    param=p.param
    if :param in fieldnames(typeof(p.model))
        param=merge(param,p.model.param)
    end
    if :kwargs in fieldnames(typeof(p.model))
        kwargs=zip(keys(p.model.kwargs),values(p.model.kwargs))
        param=merge(param,(kwargs=kwargs,))
    end

    initial = Symbol(typeof(p.initial))
    param   = Dict(pairs(param))
    solver  = Symbol(typeof(p.solver))
    data    = p.data


    ProblemSave(model,initial,param,solver,data)

end

function convert(::Type{Problem}, p :: ProblemSave)

    @show param = (;p.param...)
    @show p.model
    @show p.initial
    @show p.solver
    @show p.data.datasize,p.data.datalength

    # Reconstructs the model
    model = getfield(ShallowWaterModels, p.model)(param)
    if :kwargs in fieldnames(typeof(model))
        model = getfield(ShallowWaterModels, p.model)(param;param.kwargs...)
    end

    # Reconstructs the initial data
    mesh=Mesh(param)
    U=first(p.data.U)
    initial = Init(mesh,U[:,1],U[:,2])

    # Reconstructs the solver (may not work with other user-defined solvers)
    if p.solver == :RK4_naive
        solver = RK4_naive()
    else
        solver = getfield(ShallowWaterModels, p.solver)(param)
    end

    # Reconstructs the structure of the problem
    pb = Problem(model, initial, param; solver = solver)

    # Writes the raw data
    pb.data = p.data

    return pb

end

import JLD.save

"""
    `save(p::Problem,name::String)`

Saves the content of an object `Problem` into the file `name.jld`.
"""

function save(p::Problem,name::String)

    filename = string(name,".jld")
    JLD.save(filename, name, convert(ProblemSave,p))

end

import JLD.load

"""
    `load(name::String)`

Loads the contents of the file `name.jld` as a problem of type `:Problem`.
"""
function load(name::String)

    convert(Problem, JLD.load(string(name,".jld"), name ))

end
