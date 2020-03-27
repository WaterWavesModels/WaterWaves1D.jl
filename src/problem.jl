export Problem

"""
    `Problem( model, initial, param; solver)`

    Builds an initial-value problem which can then be solved (integrated in time) through `solve!( problem )`

# Arguments
- `model   :: AbstractModel`,  the system of equation solved. May be built, e.g., by `WhithamGreenNaghdi(param)`;
- `initial :: InitialData`, the initial data. May be buit, e.g., by `Init(η,v)` where `η` is the surface deformation and `v` the derivative of the trace of the velocity potential at the surface
- `param   :: NamedTuple`, must contain values for
    - `N`, the number of collocation points of the spatial grid
    - `L`, the half-length of the spatial grid
    - `T`, the final time of integration
    - `dt`, the timestep
    - `nr` (optional, default = `T/dt`) the number of stored data
- `solver`  :: TimeSolver (optional, default = explicit Runge-Kutta fourth order solver), the solver for time integration. May be built, e.g., by `RK4(param)` or `RK4_naive()`


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
                     param   :: NamedTuple;
                     solver = RK4(param,model)  :: TimeSolver)

        if in(:ns,keys(param))
            times = Times(param.dt, param.T; ns = param.ns)
        else
            times = Times(param.dt, param.T)
        end

        mesh  = Mesh(param)

        data  = Data(mapto(model, initial))

        new(model, initial, param, solver, times, mesh, data)

    end

end

export solve!

"""
    `solve!( problem)`

    Solves (i.e. integrates in time) an initial-value problem

The argument `problem` should be of type `:: Problem`.
It may be buit, e.g., by `Problem( model, initial, param)`

"""
function solve!(problem :: Problem)

    @show problem.param

    U = copy(last(problem.data.U))

    dt = problem.times.dt

    J = range(problem.times.ns ,stop = problem.times.Nc-1, step = problem.times.ns)

    if problem.times.ns == 1
        @showprogress 1 for j in J
            step!(problem.solver, problem.model, U, dt)
            push!(problem.data.U,copy(U))
        end

    else
        L = 1:problem.times.ns
        @showprogress 1 for j in J
            @showprogress 1 for l in L
                step!(problem.solver, problem.model, U, dt)
            end
            println()
            push!(problem.data.U,copy(U))
        end
    end

    println()

end
