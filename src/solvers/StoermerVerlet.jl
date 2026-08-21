export StoermerVerlet
export step!

@doc raw"""
$(TYPEDEF)

Störmer–Verlet solver [HairerLubichWanner2003](@citet) for canonical Hamiltonian equations.
Combination of the composition of the two symplectic Euler methods (with different equations solved implicitly).
The implicit problems are solved using explicit fixed-point iterations.

Construct an object of type `TimeSolver` to be used in `Problem(model, initial, param; solver::TimeSolver)`

Arguments can be either
0. an object of type `AbstractModel`;
1. an `Array` of size `(N,2)` where `N` is the number of collocation points;
2. a `NamedTuple` containing a key `N`.

The keyword argument `Niter` (optional, defaut value = 10) determines the number of steps in the Neumann iteration solver of the implicit step.
The keyword argument `implicit` (optional, defaut value = 1) determines which equation is first solved implicitly (must be `1` or `2`).

The function

```julia
step!(solver :: StoermerVerlet, model :: AbstractModel , U, δt)
```

performs the integration step of the symplectic Euler solver applied to solutions to the equation ``(u₁,u₂)'=(f₁,f₂)(u₁,u₂)``.

It replaces the argument ``U≈(u₁,u₂)(tₙ)`` with the next element of the recursive scheme approximating ``(u₁,u₂)(tₙ+δt)`` through the formula

```math
 \left\{\begin{array}{l}
u₁(tₙ+δt/2)≈ u₁(tₙ) + δt/2 f₁( u₁(tₙ+δt/2), u₂(tₙ) )\\
u₂(tₙ+δt/2)≈ u₂(tₙ) + δt/2 f₂( u₁(tₙ+δt/2), u₂(tₙ) )\\
u₂(tₙ+δt)≈ u₂(tₙ+δt/2) + δt/2 f₂( u₁(tₙ+δt/2), u₂(tₙ+δt) )\\
u₁(tₙ+δt)≈ u₁(tₙ+δt/2) + δt/2 f₁( u₁(tₙ+δt/2), u₂(tₙ+δt) )
  \end{array}\right.
```
(if the first equation is solved implicitly first, as defined by `implicit`).

"""
struct StoermerVerlet{T} <: TimeSolver

    U1::Vector{T}
    U2::Vector{T}
    Niter::Int
    implicit::Int
    label::String
    info::String

    function StoermerVerlet(U::Vector{Vector{T}}; Niter = 10, implicit = 1) where T
        U1 = deepcopy(U[1])
        U2 = deepcopy(U[2])
        if implicit != 1 && implicit != 2
            @warn "the keyword `implicit` must be 1 or 2. solve! will not work."
        end
        info = "Symplectic Euler time solver: equation $implicit is solved first via \
        the implicit Euler step (using the Neumann expansion with $Niter iterations) \
        and then equation $(3 - implicit) is solved via the explicit Euler step."
        label = "symplectic Euler"
        return new{T}(U1, U2, Niter, implicit, label, info)
    end

end

function StoermerVerlet(model::AbstractModel; Niter = 10, implicit = 1)
    U = model.mapto(Init(x -> 0 * x, x -> 0 * x))
    return StoermerVerlet(U; Niter = Niter, implicit = implicit)
end

function StoermerVerlet(param::NamedTuple; Niter = 10, implicit = 1)
    return StoermerVerlet([zeros(Complex{Float64}, param.N), zeros(Complex{Float64}, param.N)]; Niter = Niter, implicit = implicit)
end

function step!(
        solver::StoermerVerlet,
        model::AbstractModel,
        U,
        dt
    )


    solver.U1 .= deepcopy(U[1])
    solver.U2 .= deepcopy(U[2])

    return if solver.implicit == 1
        for i in 1:solver.Niter
            model.f1!(solver.U1, solver.U2)
            solver.U1 .= U[1] + dt / 2 * solver.U1
        end
        U[1] .= solver.U1
        model.f2!(solver.U1, solver.U2)
        U[2] .+= dt / 2 * solver.U2
        for i in 1:solver.Niter
            model.f2!(solver.U1, solver.U2)
            solver.U2 .= U[2] + dt / 2 * solver.U2
        end
        U[2] .= solver.U2
        model.f1!(solver.U1, solver.U2)
        U[1] .+= dt / 2 * solver.U1
    elseif solver.implicit == 2
        for i in 1:solver.Niter
            model.f2!(solver.U1, solver.U2)
            solver.U2 .= U[2] + dt / 2 * solver.U2
        end
        U[2] .= solver.U2
        model.f1!(solver.U1, solver.U2)
        U[1] .+= dt / 2 * solver.U1
        for i in 1:solver.Niter
            model.f1!(solver.U1, solver.U2)
            solver.U1 .= U[1] + dt / 2 * solver.U1
        end
        U[1] .= solver.U1
        model.f2!(solver.U1, solver.U2)
        U[2] .+= dt / 2 * solver.U2
    else
        error("when defining `StoermerVerlet`, the keyword `implicit` must be either 1 or 2.")
    end
end
