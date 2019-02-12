export AbstractModel
export TimeSolver
export InitialData
export Data

export Times
export Mesh
export Problem

abstract type AbstractModel end
abstract type TimeSolver end
abstract type InitialData end

mutable struct Data
    U :: Array{Array{Complex{Float64},2}}
    datasize :: Int
    datalength:: Int

    function Data( v )
        (datalength , datasize ) = size(v)
        U = [v]
        new(U, datasize, datalength)
    end
end


struct Times
    Nt   :: Int
    Nr   :: Int
    nr   :: Int
    tfin :: Float64
    dt   :: Float64
    t    :: Vector{Float64}
    tr   :: Vector{Float64}

    function Times( dt, tfin)
        t = range(0, stop=tfin, step = dt)
        Nt = length(t)
        Nr = Nt
        nr = 1
        tr = t
        new( Nt, Nr, nr, tfin, dt, t, tr)
    end

    function Times( dt, tfin, nr)
        t = range(0, stop=tfin, step = dt)
        Nt = length(t)
        tr = t[range(1, stop=Nt, step = nr)]
        Nr = length(tr)
        new( Nt, Nr, nr, tfin, dt, t, tr)
    end
end


struct Mesh

    N    :: Int64
    xmin :: Float64
    xmax :: Float64
    dx   :: Float64
    x    :: Vector{Float64}
    kmin :: Float64
    kmax :: Float64
    dk   :: Float64
    k    :: Vector{Float64}

    function Mesh( xmin :: Float64, xmax :: Float64, N :: Int64)

        dx   = (xmax-xmin)/N
        x    = zeros(Float64, N)
        x   .= range(xmin, stop=xmax, length=N+1)[1:end-1] 
        dk   = 2π/(N*dx)
        kmin = -N/2*dk
        kmax = (N/2-1)*dk
        k    = zeros(Float64, N)
        k   .= dk .* vcat(0:N÷2-1, -N÷2:-1)

        new( N, xmin, xmax, dx, x, kmin, kmax, dk, k)

    end

    function Mesh(param :: NamedTuple)

        xmin = - Float64(param.L)
        xmax =   Float64(param.L)
        N    =   param.N
        
        Mesh( xmin, xmax, N)

    end
end


"""
    Problem( model, initial, param, solver)

- model   : CGBSW or Matsuno
- initial : BellCurve
- param   : must contain N, L, T, dt for Mesh and Times, may contain additional data for Models (ϵ)
- solver  : RK4 (optional)

"""
mutable struct Problem

    model   :: AbstractModel
    initial :: InitialData
    param   :: NamedTuple
    solver  :: TimeSolver
    times   :: Times
    mesh    :: Mesh
    data    :: Data

    function Problem(model   :: AbstractModel,
                     initial :: InitialData,
                     param   :: NamedTuple,
                     solver  :: TimeSolver)

        if in(:nr,keys(param))
            times = Times(param.dt, param.T, param.nr)
        else
            times = Times(param.dt, param.T)
        end

        mesh  = Mesh(param)
        data  = Data(mapto(model,initial))

        new(model,initial,param,solver,times,mesh,data)

    end

    function Problem( model   :: AbstractModel,
                      initial :: InitialData,
                      param   :: NamedTuple)

        if in(:nr,keys(param))
            times = Times(param.dt, param.T, param.nr)
        else
            times = Times(param.dt, param.T)
        end
        mesh   = Mesh(param)
        data   = Data(mapto(model,initial))
        solver = RK4(param,model)

        new(model,initial,param,solver,times,mesh,data)

    end

end

struct ProblemSave

    model   :: AbstractModel
    initial :: InitialData
    param   :: Dict{Symbol,Any}
    solver  :: TimeSolver
    times   :: Times
    mesh    :: Mesh
    data    :: Data

    function ProblemSave(p   :: Problem)

        model   = p.model
        initial = p.initial
        param   = Dict(pairs(p.param))
        solver  = p.solver
        times   = p.times
        mesh    = p.mesh
        data    = p.data

        new(model,initial,param,solver,times,mesh,data)

    end

end

struct Problem

    model   :: AbstractModel
    initial :: InitialData
    param   :: NamedTuple
    solver  :: TimeSolver
    times   :: Times
    mesh    :: Mesh
    data    :: Data

    function Problem(p :: ProblemSave)

        dictkeys(d::Dict) = (collect(keys(d))...,)
        dictvalues(d::Dict) = (collect(values(d))...,)
        namedtuple(d::Dict{Symbol,T}) where {T} =
        NamedTuple{dictkeys(d)}(dictvalues(d))
        model   = p.model
        initial = p.initial
        param   = namedtuple(p.param)
        solver  = p.solver
        times   = p.times
        mesh    = p.mesh
        data    = p.data

        new(model,initial,param,solver,times,mesh,data)

    end

    function Problem(p :: ProblemSave, param :: NamedTuple)

        model   = p.model
        initial = p.initial
        solver  = p.solver
        times   = p.times
        mesh    = p.mesh
        data    = p.data

        new(model,initial,param,solver,times,mesh,data)
    end

end
