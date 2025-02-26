export InitialData

"""
Abstract type defining initial data from which initial-value problems can be built,
through `Problem( model :: AbstractModel, initial :: InitialData, param :: NamedTuple )`
"""
abstract type InitialData end

show(io::IO, i::InitialData) =
    try print(io,i.info)
    catch
        print(io,"Initial data: $(i.label)")
    end

function Base.:(==)(i1::InitialData, i2::InitialData)
    x=rand(2)
    i1.η(x) == i2.η(x) && i1.v(x) == i2.v(x) && i1.label == i2.label
end

export Init
"""
    Init(data ; fast, label)

Generate an initial data to be used in the function `Problem`.

`data` should contain either

- a function `η` and a function `v` (in this order);
- a `mesh` (of type `Mesh`) and two vectors representing `η(mesh.x)` and `v(mesh.x)` (in this order);
- an array of collocation points and two vectors representing `η(x)` and `v(x)` (in this order).

In the last two cases, an optional keyword argument `fast` can be set to `true`,
(default is `false`), in which case the algorithm is faster and uses less allocations, but is less precise.

In the last case, the collocation points must be regularly spaced, otherwise an `ErrorException` is raised.

If the keyword `label::String` (used to display information to the output stream)
is not provided, then it is set to the `"user-defined"`.
"""
struct Init <: InitialData

    η :: Any
    v :: Any
    label :: String

    function Init(η , v;  label = "user-defined")
        new( x->η(x) , x->v(x), label )
    end

    function Init( mesh :: Mesh, η , v ; fast = false,  label = "user-defined")
        new( x->interpolate(mesh,η,x;fast=fast) , x->interpolate(mesh,v,x;fast=fast), label )
    end

    function Init( x :: Array, η , v ; fast = false, label = "user-defined")
        if !(x[2:end].-x[2]≈x[1:end-1].-x[1])
            @error("Collocation points must be equally spaced.")
        end
        mesh=Mesh(x)
        Init(mesh, η , v ; fast = fast, label = label)
    end

end

export Init2D
"""
    Init2D(data ; fast, label)

Generate a two-dimensional initial data to be used in the function `Problem`.

`data` should contain either

- a function `η` and a function `v` (in this order);
- a `mesh` (of type `Mesh`) and two vectors representing `η(mesh.x)` and `v(mesh.x)` (in this order);
- an array of collocation points and two vectors representing `η(x)` and `v(x)` (in this order).

In the last two cases, an optional keyword argument `fast` can be set to `true`,
(default is `false`), in which case the algorithm is faster and uses less allocations, but is less precise.

In the last case, the collocation points must be regularly spaced, otherwise an `ErrorException` is raised.

If the keyword `label::String` (used to display information to the output stream)
is not provided, then it is set to the `"user-defined"`.
"""
struct Init2D <: InitialData

    η :: Any
    vx :: Any
    vy :: Any
    label :: String

    function Init2D( η , vx , vy ;  label = "user-defined")
        new( (x,y)->η(x,y), (x,y)->vx(x,y), (x,y)->vy(x,y), label )
    end

    # function Init( mesh :: Mesh, η , vx, vy ; fast = false,  label = "user-defined")
    #     new( x->interpolate(mesh,η,x;fast=fast) , x->interpolate(mesh,v,x;fast=fast), label )
    # end

    # function Init( x :: Array, η , v ; fast = false, label = "user-defined")
    #     if !(x[2:end].-x[2]≈x[1:end-1].-x[1])
    #         @error("Collocation points must be equally spaced.")
    #     end
    #     mesh=Mesh(x)
    #     Init(mesh, η , v ; fast = fast, label = label)
    # end

end
