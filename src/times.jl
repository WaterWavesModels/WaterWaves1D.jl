export Times

struct Times

    Nt   :: Int
    Nr   :: Int
    nr   :: Int
    tfin :: Float64
    dt   :: Float64
    t    :: Vector{Float64}
    tr   :: Vector{Float64}

    function Times( dt, tfin)
        t = range(0, stop=tfin, step = dt)
        Nt = length(t)
        Nr = Nt
        nr = 1
        tr = t
        new( Nt, Nr, nr, tfin, dt, t, tr)
    end

    function Times( dt, tfin, nr)
        t = range(0, stop=tfin, step = dt)
        Nt = length(t)
        tr = t[range(1, stop=Nt, step = nr)]
        Nr = length(tr)
        new( Nt, Nr, nr, tfin, dt, t, tr)
    end

end

