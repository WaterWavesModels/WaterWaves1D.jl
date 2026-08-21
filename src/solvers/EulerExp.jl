export EulerExp, EulerExp_naive
export step!

@doc raw"""
Exponential Euler solver [HochbruckOstermann2010](@citet).

Construct an object of type `TimeSolver` to be used in `Problem(model, initial, param; solver::TimeSolver)`

Arguments can be either
0. an object of type `AbstractModel`, which contains the necessary;
1. an `Array` of size `(N,datasize)` where `N` is the number of collocation points and `datasize` the number of equations solved;
2. `(param,datasize)` where `param` is a `NamedTuple` containing a key `N`, and `datasize` a integer (optional, by default `datasize=2`).

`T`: element type of the solution array (e.g. `Float64`, `ComplexF64`). `T` is set at the initialisation using the type of argument `U`.

The function

```julia
    step!(solver :: EulerExp, model :: AbstractModel , U, δt)
```

performs the integration step of the exponential Euler solver applied to solutions to the equation ``u'=D u + g(u)``.

It replaces the argument ``U≈u(tₙ)`` with the next element of the recursive scheme approximating ``u(tₙ+δt)`` through the formula

```math
u(tₙ+δt)≈e^{δt D} u(tₙ) + δt \frac{e^{δt D} - 1}{δt D} g( u(tₙ) )
```

The matrix `D` should be *diagonal* and the vector of its diagonal values provided together with the nonlinear function `g` by `model`. 
"""
struct EulerExp{T} <: TimeSolver

    U1::Vector{Vector{T}}
    D::Vector{Vector{T}}
    φ::Function
    label::String

    function EulerExp(U::Vector{Vector{T}}) where T
        U1 = deepcopy(U)
        D = deepcopy(U)
        φ(z) = (exp(z + eps()) - 1) / (z + eps())
        return new{T}(U1, D, φ, "exponential Euler")
    end

end

function EulerExp(model::AbstractModel)
    U = model.mapto(Init(x -> 0 * x, x -> 0 * x))
    return EulerExp(U)
end

function EulerExp(param::NamedTuple, systemsize = 2::Int)
    return EulerExp([Array{ComplexF64}(undef, param.N) for _ in 1:systemsize])
end

function EulerExp(datasize, systemsize = 2::Int)
    return EulerExp([Array{ComplexF64}(undef, datasize) for _ in 1:systemsize])
end


function step!(
        solver::EulerExp,
        model::AbstractModel,
        U,
        dt
    )


    [u1 .= u for (u1, u) in zip(solver.U1, U)]
    model.g!(solver.U1)
    [d .= md for (d, md) in zip(solver.D, model.D)]
    [d .*= dt for d in solver.D]
    [u .*= exp.(d) for (d, u) in zip(solver.D, U)]
    [u1 .*= solver.φ.(d) for (d, u1) in zip(solver.D, solver.U1)]
    [u1 .*= dt for u1 in solver.U1]
    return [u .+= u1 for (u, u1) in zip(U, solver.U1)]

end

"""
    EulerExp_naive()

Exponential Euler solver.

A naive version of `EulerExp`, without argument since no pre-allocation is performed.

"""
struct EulerExp_naive <: TimeSolver
    φ::Function
    label::String

    function EulerExp_naive()
        φ(z) = (exp(z + eps()) - 1) / (z + eps())
        return new(φ, "exponential Euler (naive)")
    end
end

function step!(
        solver::EulerExp_naive,
        model::AbstractModel,
        U,
        dt
    )


    U0 = deepcopy(U)
    model.g!(U0)
    [u .*= exp.(dt * d) for (d, u) in zip(model.D, U)]
    [u0 .*= solver.φ.(dt * d) for (d, u0) in zip(model.D, U0)]
    return [u .+= dt * u0 for (u, u0) in zip(U, U0)]

end
