using JLD

export ProblemSave, save, load

struct ProblemSave

    model   :: Symbol
    initial :: Symbol
    param   :: Dict
    solver  :: Symbol
    data    :: Data


end

import Base.convert

function convert(::Type{ProblemSave}, p :: Problem )

    model   = Symbol(typeof(p.model))
    initial = Symbol(typeof(p.initial))
    param   = Dict(pairs(p.param))
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

    if p.model == :CGBSW
        model = CGBSW(param)
    elseif p.model == :CGBSW_naive
        model = CGBSW_naive(param)
    elseif p.model == :Matsuno
        model = Matsuno(param)
    elseif p.model == :Matsuno_naive
        model = Matsuno_naive(param)
    elseif p.model == :Matsuno_mod_naive
        model = Matsuno_mod_naive(param)
    end

    if p.initial == :BellCurve
        initial = BellCurve(param)
    elseif p.initial == :HighFreq
        initial = HighFreq(param)
    elseif p.initial == :Random
        initial = Random(param)
    end

    if p.solver == :RK4
        solver = RK4(param)
    end

    pb = Problem(model, initial, param, solver)
    pb.data = p.data

    return pb

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
