export Midpoint
export step!

@doc raw"""
    Midpoint{T} <: TimeSolver

Implicit midpoint solver (consider `RK2` for the explicit midpoint solver).
The implicit problem is solved via Neumann iteration

Construct an object of type `TimeSolver` to be used in `Problem(model, initial, param; solver::TimeSolver)`

Arguments can be either
0. an object of type `AbstractModel`;
1. an `Array` of size `(N,2)` where `N` is the number of collocation points;
2. a `NamedTuple` containing a key `N`.

The keyword argument `Niter` (optional, defaut value = 10) determines the number of steps in the Neumann iteration solver of the implicit step.
The keyword argument `realdata` is optional, and determines whether pre-allocated vectors are real- or complex-valued.
By default, they are either determined by the model or the type of the array in case `0.` and `1.`, complex-valued in case `2.`.

The function
    `step!(solver :: Midpoint, model :: AbstractModel , U, δt)`

performs the integration step of the implicit midpoint solver applied to solutions to the equation ``u'=f(u)``.

It replaces the argument ``U≈u(tₙ)`` with the next element of the recursive scheme approximating ``u(tₙ+δt)`` through the formula

```math
u(tₙ+δt)≈ u(tₙ) + δt  f( (u(tₙ)+u(tₙ+δt))/2 )
```
or, equivalently,
```math
 \left\{\begin{array}{l}
u(tₙ+δt/2)≈ u(tₙ) + δt/2 f( u(tₙ+δt/2) ) \\
u(tₙ+δt)≈ 2 u(tₙ+δt/2) - u(tₙ) 
  \end{array}\right.
```
"""
struct Midpoint{T,N} <: TimeSolver

    U0::Vector{Array{T, N}}
    Uhalf::Vector{Array{T, N}}
    Niter::Int
    label::String


    function Midpoint(U::Vector{Array{T, N}};Niter=10) where {T, N}
        U0 = deepcopy(U)
        Uhalf = deepcopy(U)

        return new{T, N}(U0, Uhalf, Niter, "implicit midpoint")
    end

end

function Midpoint(model::AbstractModel;Niter=10)
    U = model.mapto(Init(x -> 0 * x, x -> 0 * x))
    return Midpoint(U;Niter=Niter)
end

function Midpoint(param::NamedTuple, systemsize = 2::Int;Niter=10)
    return Midpoint([zeros(ComplexF64, param.N) for _ in 1:systemsize];Niter=Niter)
end

function Midpoint(datasize, systemsize = 2::Int;Niter=10)
    return Midpoint([zeros(ComplexF64, datasize) for _ in 1:systemsize];Niter=Niter)
end

function step!(
        solver::Midpoint,
        model::AbstractModel,
        U,
        dt
    )

    solver.U0 .= deepcopy(U)
    solver.Uhalf .= deepcopy(U)

    for i in 1:solver.Niter
        model.f!(solver.Uhalf)
        solver.Uhalf .= solver.U0 + dt/2 * solver.Uhalf
    end
    U .= 2*solver.Uhalf - solver.U0
    return
end
