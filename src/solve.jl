export solve!

function solve!(problem :: Problem)

    @show problem.param

    U = similar(problem.solver.Uhat)
    U .= problem.data.U[end]

    dt = problem.times.dt

    prog = Progress(problem.times.Nt,1)

    J = range(problem.times.nr ,stop = problem.times.Nt-1, step = problem.times.nr)
    L = 1:problem.times.nr

    for j in J
        for l in L

            step!(problem.solver, problem.model, U, dt)

            next!(prog)

        end

        push!(problem.data.U,copy(U))


    end

    print("\n")

end
