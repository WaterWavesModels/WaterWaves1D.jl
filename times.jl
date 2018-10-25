struct Times
    nt   :: Int
    tfin :: Float64
    dt   :: Float64
    t    :: Vector{Float64}
    
    function Times( dt, tfin)
        t = range(0, stop=tfin, step = dt)
        nt = length(t)
        new( nt, tfin, dt, t)
    end
end
