import JLD

export ProblemLight, ProblemSave, save, load

struct ProblemLight

    model   :: String
    param   :: NamedTuple

end

struct ProblemSave

    model   :: String
    param   :: Dict

end

import Base.convert

function convert(::Type{ProblemSave}, p :: ProblemLight )

    model   = String(Symbol(typeof(p.model)))
    param   = Dict(pairs(p.param))

    ProblemSave(model,param)

end

function convert(::Type{ProblemLight}, p :: ProblemSave)

    model   = p.model
    param   = (;p.param...)

    ProblemLight(model,param)

end

import JLD.save

function save(p::ProblemLight,name::String)

    filename = string(name,".jld")
    JLD.save(filename, name, convert(ProblemSave,p))

end

import JLD.load

function load(name::String)

    convert(ProblemLight, JLD.load(string(name,".jld"), name ))

end
