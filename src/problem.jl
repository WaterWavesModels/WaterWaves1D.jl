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
	mesh    :: Mesh
    data    :: Vector{Tuple{Vector{Complex{Float64}},
			    Vector{Complex{Float64}}}}

    function Problem(model   :: AbstractModel,
         	     initial :: InitialData,
         	     param   :: Parameters,
         	     solver  :: TimeSolver)

         times = Times(param.dt, param.T)
		 mesh  = Mesh(-param.L, param.L, param.N)
         data  = []

         new(model,initial,param,solver,times,mesh,data)

    end
end

export solve!

function solve!(problem :: Problem)


    h = construct(problem.model,problem.initial)[1]
    u = construct(problem.model,problem.initial)[2]

    prog = Progress(problem.times.Nt,1)

    push!(problem.data,(h,u))

    for l in range(1,problem.times.Nt-1)

        dt = problem.times.t[l+1]-problem.times.t[l]

        step!(problem.solver, problem.model, h, u, dt)
		# TO DO : faire que (h,u) soit sol, dans un AbstractType Solution, dont le type puisse changer de modele en modele

        push!(problem.data,(h,u))
		# TO DO : raffiner times de facon a ne stocker qu'un certain nombre parmi les temps calculés.

        next!(prog)

    end

end
