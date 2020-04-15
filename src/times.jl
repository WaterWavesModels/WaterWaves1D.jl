export Times
"""
    `Times(dt, tfin; ns = 1)`

Constructs a mesh of times.

# Arguments
- `dt`  : timestep;
- `tfin`: final time;
- `ns`  : data are stored every `nr` computations (optional, default = 1).


# Return values
`t=Times(args)` is of parametric type and offers
- `t.Nc`: number of computed times;
- `t.Ns`: number of recorded (stored) times;
- `t.ns`: data are stored every `ns` computations;
- `t.tfin`: the final time;
- `t.dt`: the timestep;
- `t.tc` : the vector of computed times;
- `t.ts`: the vector of recorded times.

"""
struct Times

    Nc   :: Int
    Ns   :: Int
    ns   :: Int
    tfin :: Float64
    dt   :: Float64
    tc    :: Vector{Float64}
    ts   :: Vector{Float64}

    function Times( dt, tfin; ns = 1)
        tc = range(0, stop=tfin, step = dt)
        Nc = length(tc)
        ts = tc[range(1, stop=Nc, step = ns)]
        Ns = length(ts)
        new( Nc, Ns, ns, tfin, dt, tc, ts)
    end

end
