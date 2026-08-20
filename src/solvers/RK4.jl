export RK4, RK4_naive
export step!

@doc raw"""
    RK4(arguments;realdata)

Explicit Runge-Kutta fourth order solver.

Construct an object of type `TimeSolver` to be used in `Problem(model, initial, param; solver::TimeSolver)`

Arguments can be either
0. an object of type `AbstractModel`, typically the model you will solve with the solver;
1. an `Array` which has the size of the objects that the solver will manipulate (typically a vector of `systemsize` elements of size `N` where `N` is the number of collocation points and `systemsize` the number of solved equations);
2. a `(datasize,systemsize)` where `datasize` is the size of scalar variables (typically `N` the number of collocation points) and `datasize` (optional, by default `systemsize=2`) the number of solved equations);
3. `(param,systemsize)` where `param` is a `NamedTuple` containing a key `N` describing the number of collocation points, and `systemsize` the number of solved equations (optional, by default `systemsize=2`).

The keyword argument `realdata` is optional, and determines whether pre-allocated vectors are real- or complex-valued.
By default, they are either determined by the model or the type of the array in case `0.` and `1.`, complex-valued in case `2.`.

The function
    `step!(solver :: RK4, model :: AbstractModel , U, δt)`

performs the integration step of the standard Runge-Kutta 4 solver applied to solutions to the equation `` u'=f(u)``.

It replaces the argument ``U≈u(tₙ)`` with the next element of the recursive scheme approximating ``u(tₙ+δt)`` through the formula

```math
u(tₙ+δt)≈ u(tₙ) + δt/6  * (u₁ + 2 u₂ + 2 u₃ + u₄ )
```
where
```math
 \left\{\begin{array}{l}
u₁ = f( u(tₙ) )\\
u₂ = f( u(tₙ) + δt/2 * f( u₁ ) )\\
u₃ = f( u(tₙ) + δt/2 * f( u₂ ) )\\
u₄ = f( u(tₙ) + δt * f( u₃ ) )\\
\end{array}\right.
```
"""
struct RK4{T, N} <: TimeSolver

    U1::Vector{Array{T, N}}
    dU::Vector{Array{T, N}}
    label::String

    function RK4(U::Vector{Array{T, N}}; realdata = false) where {T, N}
        U1 = deepcopy(U)
        dU = deepcopy(U)
        if realdata
            U1 = real.(U1)
            dU = real.(dU)
            return new{Float64, N}(U1, dU, "RK4")
        else
            U1 = complex.(U1)
            dU = complex.(dU)
            return new{ComplexF64, N}(U1, dU, "RK4")
        end
    end

    function RK4(model::AbstractModel; realdata = false)
        U = model.mapto(Init(x -> 0 * x, x -> 0 * x))
        return RK4(U; realdata = realdata)
    end
    function RK4(param::NamedTuple, systemsize = 2::Int; realdata = false)
        return RK4([zeros(ComplexF64, param.N) for _ in 1:systemsize]; realdata = realdata)
    end
    function RK4(datasize, systemsize = 2::Int; realdata = false)
        return RK4([zeros(ComplexF64, datasize) for _ in 1:systemsize]; realdata = realdata)
    end
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

function step!(s::RK4, m::AbstractModel, U, dt)

    for (u1, u) in zip(s.U1, U)
        copy!(u1, u)
    end

    # k1 = f(U)
    m.f!(s.U1)

    for (du, u1) in zip(s.dU, s.U1)
        copy!(du, u1)
    end

    # k2 = f(U + dt/2*k1)
    _predict!(s.U1, U, dt / 2)
    m.f!(s.U1)
    _accumulate!(s.dU, s.U1, 2)

    # k3 = f(U + dt/2*k2)
    _predict!(s.U1, U, dt / 2)
    m.f!(s.U1)
    _accumulate!(s.dU, s.U1, 2)

    # k4 = f(U + dt*k3)
    _predict!(s.U1, U, dt)
    m.f!(s.U1)
    _accumulate!(s.dU, s.U1, 1)

    # U += dt/6 * (k1 + 2k2 + 2k3 + k4)
    return _accumulate!(U, s.dU, dt / 6)

end

"""
    RK4_naive()

Runge-Kutta fourth order solver.

A naive version of `RK4`, without argument since no pre-allocation is performed.

"""
struct RK4_naive <: TimeSolver

    label::String

    function RK4_naive()
        return new("RK4 (naive)")
    end
end

function step!(
        s::RK4_naive,
        m::AbstractModel,
        U,
        dt
    )


    U0 = deepcopy(U)
    m.f!(U0)
    U1 = deepcopy(U0)

    [u0 .= u .+ dt / 2 .* u1 for (u0, u, u1) in zip(U0, U, U1)]
    m.f!(U0)
    U2 = deepcopy(U0)

    [u0 .= u .+ dt / 2 .* u2 for (u0, u, u2) in zip(U0, U, U2)]
    m.f!(U0)
    U3 = deepcopy(U0)

    [u0 .= u .+ dt .* u3 for (u0, u, u3) in zip(U0, U, U3)]
    m.f!(U0)
    U4 = deepcopy(U0)

    return [u .+= dt / 6 .* (u1 + 2 * u2 + 2 * u3 + u4) for (u, u1, u2, u3, u4) in zip(U, U1, U2, U3, U4)]

end
