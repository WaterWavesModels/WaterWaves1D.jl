using JLD

export ProblemSave, save, load

struct ProblemSave

    model   :: Symbol
    initial :: Symbol
    param   :: Dict
    solver  :: Symbol

end

import Base.convert

function convert(::Type{ProblemSave}, p :: Problem )

    model   = Symbol(typeof(p.model))
    initial = Symbol(typeof(p.initial))
    param   = Dict(pairs(p.param))
    solver  = Symbol(typeof(p.solver))

    ProblemSave(model,initial,param, solver)

end

function convert(::Type{Problem}, p :: ProblemSave)

    @show param = (;p.param...)
    @show p.model
    @show p.initial
    @show p.solver

    if p.model == :CGBSW
        model = CGBSW(param)
    elseif p.model == :Matsuno
        model = Matsuno(param)
    end

    if p.initial == :BellCurve
        initial = BellCurve(param, 1)
    end

    if p.solver == :RK4
        solver = RK4(param)
    end

    Problem(model, initial, param, solver)

end

import JLD.save

function save(p::Problem,name::String)

    filename = string(name,".jld")
    JLD.save(filename, name, convert(ProblemSave,p))

end

import JLD.load

function load(name::String)

    convert(Problem, JLD.load(string(name,".jld"), name ))

end
