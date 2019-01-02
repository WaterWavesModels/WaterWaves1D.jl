export Times

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
