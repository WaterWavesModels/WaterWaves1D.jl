export Mesh

"""
    `Mesh(args)`

Constructs a mesh of collocation points and associated Fourier modes.

# Arguments
Can be either
- `xmin :: Float64, xmax :: Float64, N :: Int64`;
- `L :: Float64, N :: Int64`, same as above with `xmin=-L` and `xmax=L`;
- `param :: NamedTuple`, param contains `N` and `L`, then same as above.

The mesh as `N` collocation points regularly spaced between `xmin` (included) and `xmax` (excluded)

# Return values
`m=Mesh(args)` is of parametric type and offers
with
- `m.N `: number of collocation points and Fourier modes;
- `m.xmin`: minimum of the mesh;
- `m.xmax`: maximum of the mesh;
- `m.dx`: distance between two collocation points;
- `m.x`: the vector of collocation points;
- `m.kmin`: minimum of Fourier modes;
- `m.kmax`: maximum of Fourier modes;
- `m.dk`: distance between two Fourier modes;
- `m.k`: the vector of Fourier modes.

"""
struct Mesh

    N    :: Int64
    xmin :: Float64
    xmax :: Float64
    dx   :: Float64
    x    :: Vector{Float64}
    kmin :: Float64
    kmax :: Float64
    dk   :: Float64
    k    :: Vector{Float64}

    function Mesh( xmin :: Float64, xmax :: Float64, N :: Int64)

        dx   = (xmax-xmin)/N
        x    = zeros(Float64, N)
        x   .= range(xmin, stop=xmax, length=N+1)[1:end-1]
        dk   = 2π/(N*dx)
        kmin = -N/2*dk
        kmax = (N/2-1)*dk
        k    = zeros(Float64, N)
        k   .= dk .* vcat(0:N÷2-1, -N÷2:-1)

        new( N, xmin, xmax, dx, x, kmin, kmax, dk, k)

    end

    function Mesh(L :: Float64, N :: Int64)

        xmin = - L
        xmax =   L
        N    =   N

        Mesh( xmin, xmax, N)

    end

    function Mesh(param :: NamedTuple)

        xmin = - Float64(param.L)
        xmax =   Float64(param.L)
        N    =   param.N

        Mesh( xmin, xmax, N)

    end
end
