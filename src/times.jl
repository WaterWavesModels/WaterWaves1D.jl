export Times
"""
    `Times(dt, tfin; ns, Ns)`

Constructs a mesh of times.

# Arguments
- `dt`  : timestep;
- `tfin`: final time;
- `ns`  : data are stored every `ns` computations (optional, default = 1).
- `Ns`  : `Ns` data (in addition to the initial datum) are stored (optional, by default `floor( tfin/dt)).
If both `Ns` and `ns` are given, `Ns` overrules `ns`.

# Return values
`t=Times(args)` is of parametric type and offers
- `t.Nc`: number of computed times;
- `t.Ns`: number of stored times;
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
        if Ns != nothing
            ind_stored = round.(Int,range(1, stop=Nc, length = minimum([Ns+1 Nc])))
        else
            ind_stored = round.(Int,range(1, stop=Nc, step = ns))
        end
        ts = tc[ind_stored]
        Ns = length(ts)
        ns = ind_stored[2:end]-ind_stored[1:end-1]
        new( Nc, Ns, ns, tfin, dt, tc, ts)
    end

end
