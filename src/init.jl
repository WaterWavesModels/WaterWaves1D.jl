export Init
"""
    Init(data ; fast=false)

Generate an initial data to be used in the function `Problem`.

`data` should contain either

- a function `η` and a function `v` (in this order)
- a Namedtuple with a function `η` and a function `v`
- a mesh and two vectors representing η(mesh.x) and v(mesh.x) (in this order)
- a mesh and a Namedtuple with a vector `η` and a vector `v` as above
- an array of collocation points and two vectors representing `η(x)` and `v(x)` (in this order)
- a mesh and a Namedtuple with a vector `η` and a vector `v` as above

In the last four cases, an optional keyword argument `fast` can be set to `true`,
(default is `false`), in which case the algorithm is faster and uses less allocations, but is less precise.

In the last two cases, the collocation points must be regularly spaced, otherwise an `ErrorException` is raised.

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

    function Init(mesh :: Mesh, η0 , v0 ; fast = false)
        new( x->interpolate(mesh,p.η,x;fast=fast) , x->interpolate(mesh,p.v,x;fast=fast) )
    end

    function Init(mesh :: Mesh, p :: NamedTuple ; fast = false)
        new( x->interpolate(mesh,p.η,x;fast=fast) , x->interpolate(mesh,p.v,x;fast=fast) )
    end

    function Init(x :: Array, η0 , v0 ; fast = false)
        y=x[2:end]-x[1:end-1];y.-=y[1];
        if maximum(abs.(y))>8*eps(maximum(x))
            @error("Collocation points must be equally spaced.")
        end
        mesh=Mesh(x)
        Init(mesh, η0 , v0 ; fast = fast)
    end

    function Init(x :: Array, p :: NamedTuple ; fast = false)
        y=x[2:end]-x[1:end-1];y.-=y[1];
        if maximum(abs.(y))>8*eps(maximum(x))
            @error("Collocation points must be equally spaced.")
        end
        mesh=Mesh(x)
        Init(mesh, p ; fast = fast)
    end


end
