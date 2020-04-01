# -*- coding: utf-8 -*-
# # Minimal example for WhithamGreenNaghdi

using ProgressMeter
using FFTW, LinearAlgebra
using LinearMaps, IterativeSolvers
using Statistics
using Plots

# +
abstract type AbstractModel end
abstract type TimeSolver end
abstract type InitialData end

include("../src/data.jl")
include("../src/times.jl")
include("../src/mesh.jl")

# +
"""
    Problem( model, initial, param; solver)

Builds an initial-value problem which can then be solved (integrated in time) through `solve!( problem )`

# Arguments
- `model   :: AbstractModel`,  the system of equation solved. May be built, e.g., by `WhithamGreenNaghdi(param)`;
- `initial :: InitialData`, the initial data. May be buit, e.g., by `Init(η,v)` where `η` is the surface deformation and `v` the derivative of the trace of the velocity potential at the surface
- `param   :: NamedTuple`, must contain values for
    - `N`, the number of collocation points of the spatial grid
    - `L`, the half-length of the spatial grid
    - `T`, the final time of integration
    - `dt`, the timestep
    - `nr` (optional, default = `T/dt`) the number of stored data
- `solver`  :: TimeSolver (optional, default = explicit Runge-Kutta fourth order solver), the solver for time integration. May be built, e.g., by `RK4(param)` or `RK4_naive()`


"""
mutable struct Problem

    model::AbstractModel
    initial::InitialData
    param::NamedTuple
    solver::TimeSolver
    times::Times
    mesh::Mesh
    data::Data

    function Problem(
        model::AbstractModel,
        initial::InitialData,
        param::NamedTuple;
        solver = RK4(param, model)::TimeSolver,
    )

        if in(:ns, keys(param))
            times = Times(param.dt, param.T; ns = param.ns)
        else
            times = Times(param.dt, param.T)
        end

        mesh = Mesh(param)

        data = Data(mapto(model, initial))

        new(model, initial, param, solver, times, mesh, data)

    end

end

# +
"""
    solve!( problem )

Solves (i.e. integrates in time) an initial-value problem

The argument `problem` should be of type `:: Problem`.
It may be buit, e.g., by `Problem(model, initial, param)`

"""
function solve!(problem::Problem)

    @show problem.param

    U = copy(last(problem.data.U))

    dt = problem.times.dt

    Jn = range(problem.times.ns, stop = problem.times.Nc - 1, step = problem.times.ns)
    J = 1:length(Jn)

    L = 1:problem.times.ns
    for j in J
        @showprogress "Step $j / $(length(J)) ..." 1 for l in L
            step!(problem.solver, problem.model, U, dt)
        end
        push!(problem.data.U, copy(U))
        println()

    end

    println()

end

# +
"""
    Init(data)
    data should contain either
    - a function η and a function v (in this order)
    - a Namedtuple with a function η and a function v
    - a mesh and two vectors Vector{Complex{Float64}} or Vector{Float64} representing η(mesh.x) and v(mesh.x) (in this order)
    - a mesh and a Namedtuple with a vector η and a vector v as above

"""
struct Init <: InitialData

    η
    v

    function Init(mesh::Mesh, η0, v0)
        hatv = fft(v0)
        hatη = fft(η0)
        k = mesh.k
        x₀ = mesh.x[1]
        η(x::Vector{Float64}) = real.(exp.(1im * (x .- x₀) * k') * hatη / length(k))
        v(x::Vector{Float64}) = real.(exp.(1im * (x .- x₀) * k') * hatv / length(k))
        new(x -> η(x), x -> v(x))
    end


end


# +
"""
    WhithamGreenNaghdi(param;kwargs)

Defines an object of type `AbstractModel` in view of solving the initial-value problem for
the modified Green-Naghdi model proposed by V. Duchêne, S. Israwi and R. Talhouk.

# Argument
`param` is of type `NamedTuple` and must contain
- dimensionless parameters `ϵ` (nonlinearity) and `μ` (dispersion);
- numerical parameters to construct the mesh of collocation points as `mesh = Mesh(param)`
E.g. `param = ( μ  = 0.1, ϵ  = 1, N  = 2^9, L  = 10*π)`

## Keywords
- `SGN :: Bool`: if `true` computes the Serre-Green-Naghdi (SGN) instead of Whitham-Green-Naghdi (WGN) system (default is `false`);
- `iterative :: Bool`: solves the elliptic problem through GMRES if `true`, LU decomposition if `false` (default is `true`);
- `precond :: Bool`: Preconditioner of GMRES is based on WGN if `true`, SGN otherwise (default is `true`);
- `gtol :: Real`: relative tolerance of the GMRES algorithm (default is `1e-14`);
- `ktol :: Real`: tolerance of the Krasny filter (default is `0`, i.e. no filtering);
- `dealias :: Int`: dealiasing with Orlicz rule `1-dealias/(dealias+2)` (default is `0`, i.e. no dealiasing).

# Return values
This generates
1. a function `WhithamGreenNaghdi` to be called in the time-integration solver;
2. a function `mapto` which from `(η,v)` of type `InitialData` provides the  data matrix on which computations are to be executed.
3. a function `mapfro` which from such data matrix returns the Tuple of real vectors `(η,v,u)`, where

    - `η` is the surface deformation;
    - `v` is the derivative of the trace of the velocity potential;
    - `u` corresponds to the layer-averaged velocity.

"""
mutable struct WhithamGreenNaghdi <: AbstractModel

    label::String
    datasize::Int
    μ::Real
    ϵ::Real
    x::Vector{Float64}
    F₀::Vector{Complex{Float64}}
    ∂ₓ::Vector{Complex{Float64}}
    Π⅔::BitArray{1}
    Id::BitArray{2}
    FFT::Array{Complex{Float64},2}
    IFFT::Array{Complex{Float64},2}
    IFFTF₀::Array{Complex{Float64},2}
    M₀::Array{Complex{Float64},2}
    h::Vector{Complex{Float64}}
    u::Vector{Complex{Float64}}
    fftv::Vector{Complex{Float64}}
    fftη::Vector{Complex{Float64}}
    fftu::Vector{Complex{Float64}}
    hdu::Vector{Complex{Float64}}
    L::Array{Complex{Float64},2}
    Precond::Diagonal{Float64,Array{Float64,1}}
    iterate::Bool
    ktol::Real
    gtol::Real


    function WhithamGreenNaghdi(
        param::NamedTuple;
        iterate = true,
        SGN = false,
        dealias = 0,
        ktol = 0,
        gtol = 1e-14,
        precond = true,
    )
        if SGN == true
            label = string("Serre-Green-Naghdi")
        else
            label = string("Whitham-Green-Naghdi")
        end
        datasize = 2
        μ = param.μ
        ϵ = param.ϵ
        mesh = Mesh(param)
        k = mesh.k
        x = mesh.x
        x₀ = mesh.x[1]

        ∂ₓ = 1im * mesh.k
        F₁ = tanh.(sqrt(μ) * abs.(k)) ./ (sqrt(μ) * abs.(k))
        F₁[1] = 1                 # Differentiation
        if SGN == true
            F₀ = sqrt(μ) * ∂ₓ
        else
            F₀ = 1im * sqrt.(3 * (1 ./ F₁ .- 1)) .* sign.(k)
        end
        if precond == true
            Precond = Diagonal(1 ./ F₁)
        else
            Precond = Diagonal(1 .+ μ * k .^ 2) #Diagonal( ones(size(k)) )
        end
        Π⅔ = abs.(mesh.k) .<= mesh.kmax * (1 - dealias / (2 + dealias)) # Dealiasing low-pass filter
        FFT = exp.(-1im * k * (x .- x₀)')
        IFFT = exp.(1im * k * (x .- x₀)') / length(x)
        M₀ = IFFT * Diagonal(F₀) * FFT
        IFFTF₀ = IFFT * Diagonal(F₀)
        Id = Diagonal(ones(size(x)))
        h = zeros(Complex{Float64}, mesh.N)
        u, fftv, fftη, fftu, hdu = (similar(h),) .* ones(5)
        L = similar(FFT)

        new(
            label,
            datasize,
            μ,
            ϵ,
            x,
            F₀,
            ∂ₓ,
            Π⅔,
            Id,
            FFT,
            IFFT,
            IFFTF₀,
            M₀,
            h,
            u,
            fftv,
            fftη,
            fftu,
            hdu,
            L,
            Precond,
            iterate,
            ktol,
            gtol,
        )
    end
end


function (m::WhithamGreenNaghdi)(U::Array{Complex{Float64},2})
    m.fftη .= U[:, 1]
    m.h .= 1 .+ m.ϵ * ifft(m.fftη)
    m.fftv .= U[:, 2]

    m.L .= m.Id .- 1 / 3 .* m.FFT * Diagonal( 1 ./ (I * m.h)) * m.M₀ * Diagonal((I * m.h).^3) * m.IFFTF₀
    m.fftu .= m.L \ m.fftv

    m.u .= ifft(m.fftu)
    m.hdu .= m.h .* ifft(m.F₀ .* m.fftu)
    
    m.h .*= m.u
    fft!(m.h)
    U[:, 1] .= -m.∂ₓ .* m.h
    ifft!(m.fftv)
    m.u .= m.u .* m.fftv .- 1 / 2 .* m.u .^ 2 .- 1 / 2 .* m.hdu .^ 2
    fft!(m.u)
    U[:, 2] .= -m.∂ₓ .* ( m.fftη .+ m.ϵ * m.Π⅔ .* m.u )
    U[abs.(U).<m.ktol] .= 0
end

"""
    mapto(WhithamGreenNaghdi, data)
`data` is of type `InitialData`, maybe constructed by `Init(...)`.

Performs a discrete Fourier transform with, possibly, dealiasing and Krasny filter.

See documentation of `WhithamGreenNaghdi` for more details.

"""
function mapto(m::WhithamGreenNaghdi, data::InitialData)

    U = [m.Π⅔ .* fft(data.η(m.x)) m.Π⅔ .* fft(data.v(m.x))]
    U[abs.(U).<m.ktol] .= 0
    return U

end

"""
    mapfro(WhithamGreenNaghdi, data)
`data` is of type `Array{Complex{Float64},2}`, e.g. `last(p.data.U)` where `p` is of type `Problem`.

Returns `(η,v,u)`, where
- `η` is the surface deformation;
- `v` is the derivative of the trace of the velocity potential;
- `u` corresponds to the layer-averaged velocity.

Inverse Fourier transform and real part, plus solving the elliptic problem for `u`.

See documentation of `WhithamGreenNaghdi` for more details.
"""
function mapfro(m::WhithamGreenNaghdi, datum::Array{Complex{Float64},2})
    m.fftη .= datum[:, 1]
    m.h .= 1 .+ m.ϵ * ifft(m.fftη)
    m.L .= m.Id .- 1 / 3 * m.FFT * Diagonal(1 ./ m.h) * m.M₀ * Diagonal(m.h .^ 3) * m.IFFTF₀

    real(ifft(datum[:, 1])), real(ifft(datum[:, 2])), real(ifft(m.L \ datum[:, 2]))
end


# +
"""
    `RK4(param,model;k)`

    Constructs an object of type `TimeSolver` to be used in `Problem(model, initial, param; solver::TimeSolver)`

- `param ::NamedTuple` should contain a value N (number of collocation points)
- `model ::AbstractModel` is optional and determines the number of equations solved
- `k=2   ::Int` is optional (default = 2) and determines the number of equations solved

"""
mutable struct RK4 <: TimeSolver

    Uhat::Array{Complex{Float64},2}
    dU::Array{Complex{Float64},2}

    function RK4(param::NamedTuple, model::AbstractModel)

        Uhat = zeros(Complex{Float64}, (param.N, model.datasize))
        dU = zeros(Complex{Float64}, (param.N, model.datasize))

        new(Uhat, dU)

    end

    function RK4(param::NamedTuple; k = 2::Int)

        Uhat = zeros(Complex{Float64}, (param.N, k))
        dU = zeros(Complex{Float64}, (param.N, k))

        new(Uhat, dU)

    end

end

# +
function step!(s::RK4, f!::AbstractModel, U::Array{Complex{Float64},2}, dt::Float64)


    s.Uhat .= U
    f!(s.Uhat)
    s.dU .= s.Uhat

    s.Uhat .= U .+ dt / 2 .* s.Uhat
    f!(s.Uhat)
    s.dU .+= 2 .* s.Uhat

    s.Uhat .= U .+ dt / 2 .* s.Uhat
    f!(s.Uhat)
    s.dU .+= 2 .* s.Uhat

    s.Uhat .= U .+ dt .* s.Uhat
    f!(s.Uhat)
    s.dU .+= s.Uhat

    U .+= dt / 6 .* s.dU

end
# -

"""
	solveandplot()

Solves an initial value problem with WhithamGreenNaghdi model.
See parameters below.

"""
function solveandplot(param)
    # Set up
    

    mesh = Mesh(param) # construct mesh of collocation points, Fourier modes, etc.
    η = exp.(-(mesh.x) .^ 2)
    v = 0 * η  # Initial data

    init = Init(mesh, η, v)
    model = WhithamGreenNaghdi(param; iterate = false) # iterate = true for GMRES iterate solver
    problem = Problem(model, init, param)

    # Compute
    @time solve!(problem)

    # Plot
    fftηfin = last(problem.data.U)[:, 1]
    plt = plot(layout = (1, 2))
    plot!(
        plt[1, 1],
        fftshift(mesh.k),
        fftshift(log10.(abs.(fftηfin)));
        title = "frequencies",
    )
    plot!(plt[1, 2], mesh.x, real.(ifft(fftηfin)); title = "physical space")
    display(plt)
end

# +
param = (
        μ = 1,
        ϵ = 1,
        N = 2^8,   # number of collocation points
        L = 10 * π,   # mesh of size 2L
        T = 0.1,      # final time
        dt = 1 / 10^4, # timestep
        ns = 1000,      # stores data every ns time steps
    )

solveandplot(param)
# -


