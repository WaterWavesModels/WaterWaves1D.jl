export Problem

"""
    Problem( model, initial, param ; solver, label)

Build an initial-value problem which can then be solved (i.e. integrated in time) through `solve!( problem )`

# Arguments
- `model   :: AbstractModel`,  the system of equation solved.
May be built, e.g., by `WaterWaves(param)`;
- `initial :: InitialData`, the initial data.
May be buit, e.g., by `Init(η,v)` where
`η` is the surface deformation and `v` the derivative of the trace of the velocity potential at the surface;
- `times :: Times` is the time grid, and may be built using the function `Times`. 
Alternatively one can simply provide a `NamedTuple` with
    - `T`, the final time of integration
    - `dt`, the timestep
    - optionally, `Ns` the number of computed data or `ns` for storing data every `ns` computation steps (by default, every computed data is stored).

## Optional keyword arguments
- `solver :: TimeSolver`, the solver for time integration (default is explicit Runge-Kutta fourth order solver).
May be built, e.g., by `RK4(model)` or `RK4_naive()`.
- `label   :: String ` is used in future references (e.g. `plot_solution`).
- Information are not printed if `verbose = false` (default is `true`).

"""
mutable struct Problem

    model   :: AbstractModel
    initial :: InitialData
    solver  :: TimeSolver
    times   :: Times
    data    :: Data
    label   :: String

    function Problem(model   :: AbstractModel,
        initial :: InitialData,
        times   :: Times;
        solver = RK4(model)  :: TimeSolver,
        label = nothing
        )
        
        if isnothing(label)   label = model.label end

        U = model.mapto(initial)
        data  = Data(U)

        new(model, initial, solver, times, data, label)

    end

    function Problem(model   :: AbstractModel,
        initial :: InitialData,
        param   :: NamedTuple;
        solver = RK4(model)  :: TimeSolver,
        label = nothing
        )


        if in(:Ns,keys(param))
        times = Times(param; Ns = param.Ns)
        elseif in(:ns,keys(param))
        times = Times(param; ns = param.ns)
        else
        times = Times(param)
        end

        p = Problem(model, initial, times; solver = solver, label = label)

        new(p.model, p.initial, p.solver, p.times, p.data, p.label)

    end

end

show(io::IO, p::Problem) =
    print(io,"Initial-value problem with...\n\
    ├─Model: $(p.model.label).\n├─Initial data: $(p.initial.label).\n├─Solver: $(p.solver.label).\n\
    └─Grid of times: $(p.times.Nc) computed times on [0, $(p.times.tfin)] (timestep dt=$(p.times.dt)), \
    among which $(p.times.Ns) will be stored.")

#     Discretized with $(p.mesh.N) collocation points on [$(p.mesh.xmin), $(p.mesh.xmax)],\n\

export solve!

"""
    solve!( problem :: Problem; verbose=true )

Solve (i.e. integrate in time) an initial-value problem

The argument `problem` should be of type `Problem`.
It may be buit, e.g., by `Problem(model, initial, param)`

Information are not printed if keyword argument `verbose = false` (default is `true`).

"""
function solve!(problem :: Problem; verbose=true::Bool)

    ci = get(ENV, "CI", nothing) == "true"

    U = deepcopy(last(problem.data.U))

    dt     = problem.times.dt
    solver = problem.solver
    model  = problem.model
    data   = problem.data.U

    if verbose == true
        @info "Now solving the initial-value problem $(problem.label)\n\
            with timestep dt=$(problem.times.dt), final time T=$(problem.times.tfin),\n\
            and N=$(length(U[1])) collocation points."
    end


    if problem.times.tc == problem.times.ts
        pbar = Progress(problem.times.Ns-1; enabled = !ci)
        for j in 1:problem.times.Ns-1
            step!(solver, model, U, dt)
            push!(data,deepcopy(U))
            next!(pbar)
        end

    elseif length(problem.times.ts) > 25
        pbar = Progress(problem.times.Ns-1; enabled = !ci)
        for j in 1:problem.times.Ns-1
            for l in 1:problem.times.ns[j]
                step!(solver, model, U, dt)
            end
            push!(data,deepcopy(U))
            next!(pbar)
        end
    else
        for j in 1:problem.times.Ns-1
            pbar = Progress(problem.times.ns[j], desc=string("Step ",j,"/",problem.times.Ns-1,"...") ; enabled = !ci)
            for l in 1:problem.times.ns[j]
                step!(solver, model, U, dt)
                next!(pbar)
            end
            push!(data,deepcopy(U))
            println()

        end
    end

    println()

end



using Base.Threads

"""
    solve!( problems; verbose=true )

Solve (i.e. integrate in time) a collection of initial-value problems.

The argument `problems` should be a collection (list, array...) of elements of type `Problem`.

Information are not printed if keyword argument `verbose = false` (default is `true`).
"""
function solve!(problems; verbose=true::Bool)

    ci = get(ENV, "CI", nothing) == "true"

    # Set up
    U=[];nsteps=0;flag = false;
    for problem in problems
        for p in problems # Check whether multi-threading is at risk, using same allocations
            if !(p===problem) && ( p.model===problem.model || p.solver===problem.solver )
                flag = true
            end
        end
        push!(U,deepcopy(last(problem.data.U)))
        nsteps+=problem.times.Ns-1
    end
    if flag
        @warn "`solve!(problems)` uses multi-threading. \
        It is ill-advised to use it with problems sharing the same model or the same solver. \
        Re-define model/solver when building your problems, or solve them independently \
        via `for problem in problems solve!(problem) end`."
    end

    pg=Progress(nsteps;dt=1,enabled = !ci)
    @threads for i in 1:length(problems)
        if verbose == true
            @info "Now solving the initial-value problem $(problems[i].label)\n\
                with timestep dt=$(problems[i].times.dt), final time T=$(problems[i].times.tfin),\n\
                and N=$(length(problems[i].data.U[1][1])) collocation points."
        end

        if problems[i].times.tc == problems[i].times.ts

            for j in 1:problems[i].times.Ns-1
                step!(problems[i].solver, problems[i].model, U[i], problems[i].times.dt)
                push!(problems[i].data.U,deepcopy(U[i]))
                next!(pg)
            end

        else
            for j in 1:problems[i].times.Ns-1
                for l in 1:problems[i].times.ns[j]
                    step!(problems[i].solver, problems[i].model, U[i], problems[i].times.dt)
                end
                push!(problems[i].data.U,deepcopy(U[i]))
                next!(pg)
            end
        end
        if verbose == true
            @info "Done solving the problem $(problems[i].label)."
        end

    end

    println()

end


Base.:(≈)(x::NamedTuple{N,T}, y::NamedTuple{N2,T2}) where {N,T,N2,T2} =
    length(N) === length(union(N,N2)) && all(k->getfield(x,k) == getfield(y,k), keys(x))


Base.:(==)(p1::Problem, p2::Problem) =
    p1.model == p2.model &&
    p1.initial == p2.initial &&
    p1.solver == p2.solver &&
    p1.times == p2.times &&
    p1.data == p2.data &&
    p1.label == p2.label

Base.:(≈)(p1::Problem, p2::Problem ; atol::Real=0, rtol::Real= √eps()) =
    all(isapprox.(p1.model.mapfro(p1.data.U[end]) , p2.model.mapfro(p2.data.U[end]),atol=atol,rtol=rtol))
