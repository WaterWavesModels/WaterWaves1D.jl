export Mesh

"""
    Mesh(args)

Constructs a mesh of collocation points and associated Fourier modes.

# Arguments
Can be either
- `xmin`, `xmax`, and `N`; or
- `L`, `N` (same as above with `xmin=-L` and `xmax=L`); or
- `param :: NamedTuple`, a `NamedTuple` containing `N` and `L` or `xmin` and `xmax`, then same as above; or
- `x` a vector of regularly spaced collocation points`.

The mesh as `N` collocation points regularly spaced between `xmin` (included) and `xmax` (excluded)

# Return values
`m=Mesh(args)` is of parametric type and offers
with
- `m.N `: number of collocation points and Fourier modes;
- `m.xmin`: minimum of the mesh (included in the vector of collocation points);
- `m.xmax`: maximum of the mesh (excluded in the vector of collocation points);
- `m.dx`: distance between two collocation points;
- `m.x`: the vector of collocation points;
- `m.kmin`: minimum of Fourier modes (included in the vector of Fourier modes);
- `m.kmax`: maximum of Fourier modes (included in the vector of Fourier modes);;
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

    function Mesh( xmin , xmax , N :: Int64)

        dx   = (xmax-xmin)/N
        x    = zeros(Float64, N)
        x   .= range(xmin, stop=xmax, length=N+1)[1:end-1]
        dk   = 2π/(N*dx)
        kmin = -N÷2*dk
        kmax = (N-1)÷2*dk
        k    = zeros(Float64, N)
        k   .= dk .* vcat(0:(N-1)÷2, -N÷2:-1)

        new( N, xmin, xmax, dx, x, kmin, kmax, dk, k)

    end

    function Mesh(L , N :: Int64)

        xmin = - L
        xmax =   L
        N    =   N

        Mesh( xmin, xmax, N)

    end

    function Mesh( x )

        N    =   length(x)
        dx   =   (x[end]-x[1])/(N-1)
        xmin = x[1]
        xmax = x[end]+dx

        Mesh( xmin, xmax, N)

    end

    function Mesh(param :: NamedTuple)

        if :xmin in keys(param) && :xmax in keys(param)
            xmin = param.xmin
            xmax = param.xmax
        elseif :L in keys(param)
            xmin = - param.L
            xmax =   param.L
        else
            @error("the NamedTuple must contain the field `L`, or `xmin` and `xmax`, when defining a mesh.")
        end
        N    =   param.N

        Mesh( xmin, xmax, N)

    end
end

show(io::IO, m::Mesh) =
    print(io,"One-dimensional grid of $(m.N) collocation points on [$(m.xmin), $(m.xmax)].\n\
    Grid spacing dx=$(m.dx).")

function dump( h5file :: String, mesh :: Mesh )

    h5write(joinpath(h5file * ".h5"), "/mesh/xmin", mesh.xmin)
    h5write(joinpath(h5file * ".h5"), "/mesh/xmax", mesh.xmax)
    h5write(joinpath(h5file * ".h5"), "/mesh/N", mesh.N)

end

function load_mesh( h5file :: String )

    xmin = h5read(joinpath(h5file * ".h5"), "/mesh/xmin")
    xmax = h5read(joinpath(h5file * ".h5"), "/mesh/xmax")
    N = h5read(joinpath(h5file * ".h5"), "/mesh/N")

    return Mesh( xmin, xmax, N)

end

Base.:(==)(a::Mesh, b::Mesh) = ( a.N == b.N && a.xmin == a.xmin && a.xmax == b.xmax )
