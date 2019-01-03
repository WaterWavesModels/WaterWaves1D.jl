export solve!

function solve!(problem :: Problem)
	@show problem.param

    U = copy(problem.data.U[end])

    prog = Progress(problem.times.Nt,1)

    for l in range(1,problem.times.Nt-1)

        dt = problem.times.dt

        step!(problem.solver, problem.model, U, dt)
		# TO DO : faire que (h,u) soit sol, dans un AbstractType Solution, dont le type puisse changer de modele en modele

        push!(problem.data.U,copy(U))
		# TO DO : raffiner times de facon a ne stocker qu'un certain nombre parmi les temps calculés.

        next!(prog)

    end

	print("\n")

end
