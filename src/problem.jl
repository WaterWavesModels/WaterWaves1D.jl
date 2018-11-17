export Problem

"""
    Problem( model, initial, param, solver)

- model   : CGBSW or Matsuno
- initial : BellCurve
- param   : Mesh, Frequency, epsilon
- solver  : RK4

"""
struct Problem

    model   :: AbstractModel
    initial :: InitialData
    param   :: Parameters
    solver  :: TimeSolver
    times   :: Times
    data    :: Vector{Tuple{Vector{Complex{Float64}},
			    Vector{Complex{Float64}}}}

    function Problem(model   :: AbstractModel,
         	     initial :: InitialData,
         	     param   :: Parameters,
         	     solver  :: TimeSolver)

         times = Times(param.dt, param.T)
         data  = []

         new(model,initial,param,solver,times,data)

    end
end

export solve!

function solve!(problem :: Problem)


    h = problem.model.Pi .* fft(problem.initial.h)
    u = problem.model.Pi .* fft(problem.initial.u)

    prog = Progress(problem.times.Nt,1)

    push!(problem.data,(h,u))

    for l in range(1,problem.times.Nt-1)

        dt = problem.times.t[l+1]-problem.times.t[l]

        step!(problem.solver, problem.model, h, u, dt)

        push!(problem.data,(h,u))

        next!(prog)

    end

end
