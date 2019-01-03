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
abstract type Data end



struct Times
    Nt   :: Int
    tfin :: Float64
    dt   :: Float64
    t    :: Vector{Float64}

    function Times( dt, tfin)
        t = range(0, stop=tfin, step = dt)
        Nt = length(t)
        new( Nt, tfin, dt, t)
    end
end


struct Mesh

    N   :: Int
    xmin :: Float64
    xmax :: Float64
    dx   :: Float64
    x    :: Vector{Float64}
    kmin :: Float64
    kmax :: Float64
    dk   :: Float64
    k    :: Vector{Float64}

    function Mesh( xmin, xmax, N)
        dx = (xmax-xmin)/N
        x = range(xmin, stop=xmax, length=N+1)[1:end-1]
        dk = 2π/(N*dx)
        kmin = -N/2*dk
        kmax = (N/2-1)*dk
        k = [range(0, length=N ÷ 2, step = dk) ; range(kmin, length=N ÷ 2, step = dk) ]
        new( N, xmin, xmax, dx, x, kmin, kmax, dk, k)
    end
end


"""
    Problem( model, initial, param, solver)

- model   : CGBSW or Matsuno
- initial : BellCurve
- param   : must contain N, L, T, dt for Mesh and Times, may contain additional data for Models (ϵ)
- solver  : RK4 (optional)

"""
struct Problem

    model   :: AbstractModel
    initial :: InitialData
    param   :: NamedTuple
    solver  :: TimeSolver
    times   :: Times
	mesh    :: Mesh
    data    :: Array{Array{Complex{Float64},2}}

    function Problem(model   :: AbstractModel,
         	     initial :: InitialData,
         	     param   :: NamedTuple,
         	     solver  :: TimeSolver)

         times = Times(param.dt, param.T)
		 mesh  = Mesh(-param.L, param.L, param.N)
         data  = [mapto(model,initial)]

         new(model,initial,param,solver,times,mesh,data)

    end

	function Problem(model   :: AbstractModel,
         	     initial :: InitialData,
         	     param   :: NamedTuple)

         times = Times(param.dt, param.T)
		 mesh  = Mesh(-param.L, param.L, param.N)
         data  = [mapto(model,initial)]
		 solver= RK4(param,2)

         new(model,initial,param,solver,times,mesh,data)

    end
end
