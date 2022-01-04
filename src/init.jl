export Init

"""
    Init(data)

Generate an initial data to be used in the function `Problem`.

`data` should contain either

- a function `η` and a function `v` (in this order)
- a Namedtuple with a function `η` and a function `v`
- a mesh and two vectors representing η(mesh.x) and v(mesh.x) (in this order)
- a mesh and a Namedtuple with a vector `η` and a vector `v` as above
- an array of collocation points and two vectors representing `η(x)` and `v(x)` (in this order)
- a mesh and a Namedtuple with a vector `η` and a vector `v` as above

"""
struct Init <: InitialData

    η
    v

    function Init(η , v)
        new( x->η(x) , x->v(x) )
    end

    function Init(p :: NamedTuple)
        new( x->p.η(x) , x->p.v(x) )
    end

    function Init(mesh :: Mesh, η0 , v0)
        hatv=fft(v0)
        hatη=fft(η0)
        k = mesh.k
        x₀ = mesh.x[1]
        function η( x :: Vector{Float64} )
            return real.(exp.(1im*(x.-x₀)*k')*hatη/length(k))
        end
        function v( x :: Vector{Float64} )
            return real.(exp.(1im*(x.-x₀)*k')*hatv/length(k))
        end

        new( x->η(x) , x->v(x) )
    end

    function Init(mesh :: Mesh, p :: NamedTuple)
        hatv=fft(p.v)
        hatη=fft(p.η)
        k = mesh.k
        x₀ = mesh.x[1]
        function η( x :: Vector{Float64} )
            return real.(exp.(1im*(x.-x₀)*k')*hatη/length(k))
        end
        function v( x :: Vector{Float64} )
            return real.(exp.(1im*(x.-x₀)*k')*hatv/length(k))
        end

        new( x->η(x) , x->v(x) )
    end

    function Init(x :: Array, η0 , v0)
        @warn("Collocation points must be equally spaced.")
        mesh=Mesh(x)
        Init(mesh, η0 , v0 )
    end

    function Init(x :: Array, p :: NamedTuple)
        @warn("Collocation points must be equally spaced.")
        mesh=Mesh(x)
        Init(mesh, p )
    end


end
