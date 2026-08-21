export Euler, Euler_naive
export step!

@doc raw"""
    Euler

Explicit Euler solver.

Construct an object of type `TimeSolver` to be used in `Problem(model, initial, param; solver::TimeSolver)`

Arguments can be either
0. an object of type `AbstractModel`;
1. an `Array` of size `(N,datasize)` where `N` is the number of collocation points and `datasize` the number of equations solved;
2. `(param,datasize)` where `param` is a `NamedTuple` containing a key `N`, and `datasize` a integer (optional, by default `datasize=2`).

The function

```julia
   step!(solver :: EulerExp, model :: AbstractModel , U, δt)
```

performs the integration step of the explicit Euler solver applied to solutions to the equation ``u'=f(u)``.

It replaces the argument ``U≈u(tₙ)`` with the next element of the recursive scheme approximating ``u(tₙ+δt)`` through the formula

```math
u(tₙ+δt)≈ u(tₙ) + δt f( u(tₙ) )
```

"""
struct Euler{T,N} <: TimeSolver

    U1::Vector{Array{T,N}}
    label::String

    function Euler(U::Vector{Array{T,N}}) where {T,N}
        U1 = deepcopy(U)
        return new{T,N}(U1, "Euler")
    end
end

function Euler(model::AbstractModel)
    U = model.mapto(Init(x -> 0 * x, x -> 0 * x))
    return Euler(U)
end

function Euler(param::NamedTuple, systemsize = 2::Int)
    return Euler([Array{Complex{Float64}}(undef, param.N) for _ in 1:systemsize])
end

function Euler(datasize, systemsize = 2::Int)
    return Euler([Array{Complex{Float64}}(undef, datasize) for _ in 1:systemsize])
end

function step!(
        solver::Euler,
        model::AbstractModel,
        U,
        dt
    )


    [u1 .= u for (u1, u) in zip(solver.U1, U)]
    model.f!(solver.U1)
    [u1 .*= dt for u1 in solver.U1]
    return [u .+= u1 for (u, u1) in zip(U, solver.U1)]

end

"""
    Euler_naive()

Explicit Euler solver.

A naive version of `Euler`, without argument since no pre-allocation is performed.

"""
struct Euler_naive <: TimeSolver
    label::String

    function Euler_naive()
        return new("Euler (naive)")
    end
end

function step!(
        s::Euler_naive,
        model::AbstractModel,
        U,
        dt
    )


    U0 = deepcopy(U)
    model.f!(U0)
    return [u .+= dt * u0 for (u, u0) in zip(U, U0)]

end
