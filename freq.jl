struct Freq
    nk   :: Int
    kmin :: Float64
    kmax :: Float64
    dk   :: Float64
    k    :: Vector{Float64}
    
    function Freq( dx, nk)
        dk = 2π/(nk*dx)
        kmin = -nk/2*dk
        kmax = (nk/2-1)*dk
        k = [range(0, length=nk ÷ 2, step = dk) ; range(kmin, length=nk ÷ 2, step = dk) ]
        new( nk, kmin, kmax, dk, k)
    end
end
