using JLD

abstract type AbstractSolver end

struct RK4 <: AbstractSolver
    param
end

struct Problem

    param   :: NamedTuple
    solver  :: AbstractSolver

end

struct ProblemSave

    param   :: Dict
    solver  :: Symbol

    ProblemSave( param, solver) = new( param, solver)

end

import Base.convert

function convert(::Type{ProblemSave}, p :: Problem )

    param   = Dict(pairs(p.param))
    solver  = Symbol(typeof(p.solver))

    ProblemSave(param, solver)

end

function convert(::Type{Problem}, p :: ProblemSave)

    @show param = (;p.param...)
    @show p.solver

    if p.solver == :RK4
        solver = RK4(param)
    end

    Problem(param, solver)

end

import JLD.save

function save(p::Problem, name::String)

    filename = string(name,".jld")
    JLD.save(filename, name, convert(ProblemSave, p))

end

import JLD.load

function load(name::String)

    convert(Problem, JLD.load(string(name,".jld"), name ))

end

param = ( ϵ  = 1/2,
          N  = 2^10,
          L  = 10,
          T  = 5,
          dt = 0.01,
          theta = 1)

problem = Problem( param, RK4(param) )

p_dump  = convert(ProblemSave, problem )
p_load  = convert(Problem, p_dump )

save(problem, "problem")

p_load = load("problem")
