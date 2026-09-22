export RK2, RK2_naive
export step!

@doc raw"""
    struct RK2{T, N} <: TimeSolver

Explicit Runge-Kutta second order solver.

Construct an object of type `TimeSolver` to be used in `Problem(model, initial, param; solver::TimeSolver)`

Arguments can be either
0. an object of type `AbstractModel`, typically the model you will solve with the solver;
1. an `Array` which has the size of the objects that the solver will manipulate (typically a vector of `systemsize` elements of size `N` where `N` is the number of collocation points and `systemsize` the number of solved equations);
2. a `(datasize,systemsize)` where `datasize` is the size of scalar variables (typically `N` the number of collocation points) and `datasize` (optional, by default `systemsize=2`) the number of solved equations);
3. `(param,systemsize)` where `param` is a `NamedTuple` containing a key `N` describing the number of collocation points, and `systemsize` the number of solved equations (optional, by default `systemsize=2`).

Optionally, the keword argument `α` may be provided, and determines coefficients of the Butcher tableau. 
The default parameter is `α=1/2`, corresponding to the midpoint method. Other standard choices are `α=1` (Heun's method) and `α=2/3` (Ralston's method). 

The function

```julia
step!(solver :: RK2, model :: AbstractModel , U, δt)
```

performs the integration step of the Runge-Kutta 2 solver applied to solutions to the equation `` u'=f(u)``.

It replaces the argument ``U≈u(tₙ)`` with the next element of the recursive scheme approximating ``u(tₙ+δt)`` through the formula

```math
u(tₙ+δt)≈ u(tₙ) + δt   ((1-1/(2α))u₁ + 1/(2α) u₂  )
```
where
```math
 \left\{\begin{array}{l}
u₁ = f( u(tₙ) )\\
u₂ = f( u(tₙ) + α δt u₁ )
\end{array}\right.
```
"""
struct RK2{T, N} <: TimeSolver

    U1::Vector{Array{T, N}}
    dU::Vector{Array{T, N}}
    label::String
    α::Float64

    function RK2(U::Vector{Array{T, N}};α=1/2) where {T, N}
        U1 = deepcopy(U)
        dU = deepcopy(U)
        return new{T, N}(U1, dU, "RK2",α)
    end

end

function RK2(model::AbstractModel;α=1/2)
    U = model.mapto(Init(x -> 0 * x, x -> 0 * x))
    return RK2(U;α=α)
end

function RK2(param::NamedTuple, systemsize = 2::Int;α=1/2)
    return RK2([zeros(ComplexF64, param.N) for _ in 1:systemsize];α=α)
end

function RK2(datasize, systemsize = 2::Int;α=1/2)
    return RK2([zeros(ComplexF64, datasize) for _ in 1:systemsize];α=α)
end

@inline function _predict!(U1, U, coef)
    for (u1, u) in zip(U1, U)
        u1 .= u .+ coef .* u1
    end
    return
end

@inline function _accumulate!(dU, U1, coef)
    for (du, u1) in zip(dU, U1)
        du .+= coef .* u1
    end
    return
end

function step!(s::RK2, m::AbstractModel, U, dt)

    for (u1, u) in zip(s.U1, U)
        copy!(u1, u)
    end

    # u₁ = f(U)
    m.f!(s.U1)

    for (du, u1) in zip(s.dU, s.U1)
        copy!(du, u1)
        du.*=1-1/(2*s.α)
    end

    # u₂ = f(U + dt*α*u₁)
    _predict!(s.U1, U, dt * s.α)
    m.f!(s.U1)
    _accumulate!(s.dU, s.U1, 1/(2*s.α))


    # U += dt/6 * ((1-1/(2α))u₁ + 1/(2α)u₂)
    return _accumulate!(U, s.dU, dt )

end

"""
    RK2_naive()

Runge-Kutta fourth order solver.

A naive version of `RK2`, without argument since no pre-allocation is performed.

"""
struct RK2_naive <: TimeSolver

    label::String
    α::Float64

    function RK2_naive(;α=1/2)
        return new("RK2_naive (naive)",α)
    end
end

function step!(
        s::RK2_naive,
        m::AbstractModel,
        U,
        dt
    )


    U0 = deepcopy(U)
    m.f!(U0)
    U1 = deepcopy(U0)

    [u0 .= u .+ dt *s.α .* u1 for (u0, u, u1) in zip(U0, U, U1)]
    m.f!(U0)
    U2 = deepcopy(U0)

    return [u .+= dt .* ( (1-1/(2*s.α)) * u1 + 1/(2*s.α) * u2 ) for (u, u1, u2) in zip(U, U1, U2)]

end
