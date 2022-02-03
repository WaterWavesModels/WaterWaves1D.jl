export Times
"""
    Times(param; ns, Ns)

Constructs a mesh of times, to be used in initial-value problems (see `Problem`).

# Arguments
`param` is either
- `dt,T` with `dt` the timestep and `T` the final time of comuptation; or
- a `NamedTuple` containing `dt` and `T`

## Optional keyword arguments
- `ns`  : data are stored every `ns` computations (optional, default = 1).
- `Ns`  : `Ns` data (in addition to the initial datum) are stored (optional, by default `floor( tfin/dt)).
If both `Ns` and `ns` are given, `Ns` overrules `ns`.

# Return values
`t=Times(args)` is of parametric type and offers
- `t.Nc`: number of computed times (including initial datum);
- `t.Ns`: number of stored times (including initial datum);
- `t.ns`: number of computed times between two stored times;
- `t.tfin`: the final time;
- `t.dt`: the timestep;
- `t.tc` : the vector of computed times;
- `t.ts`: the vector of stored times.

"""
struct Times

    Nc   :: Int
    Ns   :: Int
    ns   :: Vector{Int}
    tfin :: Float64
    dt   :: Float64
    tc    :: Vector{Float64}
    ts   :: Vector{Float64}

    function Times( dt, tfin; ns = 1, Ns=nothing)
        tc = range(0, stop=tfin, step = dt)
        Nc = length(tc)
        if !isnothing(Ns)
            ind_stored = round.(Int,range(1, stop=Nc, length = minimum([Ns+1 Nc])))
        else
            ind_stored = round.(Int,range(1, stop=Nc, step = ns))
        end
        ts = tc[ind_stored]
        Ns = length(ts)
        ns = ind_stored[2:end]-ind_stored[1:end-1]
        new( Nc, Ns, ns, tfin, dt, tc, ts)
    end

    function Times( param :: NamedTuple ; ns = 1, Ns=nothing)
        Times(param.dt,param.T;ns=ns,Ns=Ns)
    end


end

show(io::IO, t::Times) =
    print(io,"Mesh of times on [0, $(t.tfin)], with timestep dt=$(t.dt).\n\
    There will be $(t.Nc) computed times (including initial data),\n\
    among which $(t.Ns) will be stored.")
