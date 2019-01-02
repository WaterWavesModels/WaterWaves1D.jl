export Mesh

struct Mesh
    
    N   :: Int
    xmin :: Float64
    xmax :: Float64
    dx   :: Float64
    x    :: Vector{Float64}
    kmin :: Float64
    kmax :: Float64
    dk   :: Float64
    k    :: Vector{Float64}
    
    function Mesh( xmin, xmax, N)
        dx = (xmax-xmin)/N
        x = range(xmin, stop=xmax, length=N+1)[1:end-1]
        dk = 2π/(N*dx)
        kmin = -N/2*dk
        kmax = (N/2-1)*dk
        k = [range(0, length=N ÷ 2, step = dk) ; range(kmin, length=N ÷ 2, step = dk) ]
        new( N, xmin, xmax, dx, x, kmin, kmax, dk, k)
    end
end
