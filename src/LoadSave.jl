export ProblemSave, save, load

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
    #
    # if p.model == :CGBSW
    #     model = CGBSW(param)
    # elseif p.model == :CGBSW_naive
    #     model = CGBSW_naive(param)
    # elseif p.model == :Matsuno
    #     model = Matsuno(param)
    # elseif p.model == :Matsuno_naive
    #     model = Matsuno_naive(param)
    # elseif p.model == :Matsuno_mod_naive
    #     model = Matsuno_mod_naive(param)
    # elseif p.model == :Boussinesq
    #     model = Boussinesq(param)
    # elseif p.model == :fdBoussinesq
    #     model = fdBoussinesq(param)
    # elseif p.model == :WaterWaves
    #     model = WaterWaves(param)
    # elseif p.model == :WhithamGreenNaghdi
    #     model = WhithamGreenNaghdi(param)
    # elseif p.model == :WhithamGreenNaghdiSym
    #     model = WhithamGreenNaghdiSym(param)
    # elseif p.model == :WhithamGreenNaghdiKlein
    #     model = WhithamGreenNaghdiKlein(param)
    # end
    model = getfield(ShallowWaterModels, p.model)(param)
    if :kwargs in fieldnames(typeof(model))
        model = getfield(ShallowWaterModels, p.model)(param;param.kwargs...)
    end

    mesh=Mesh(param)
    U=first(p.data.U)
    initial = Init(mesh,U[:,1],U[:,2])
    #
    # if p.solver == :RK4
    #     solver = RK4(param)
    # end
    solver = getfield(ShallowWaterModels, p.solver)(param)


    pb = Problem(model, initial, param; solver = solver)
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
