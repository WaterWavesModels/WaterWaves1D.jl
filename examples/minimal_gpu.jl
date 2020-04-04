# -*- coding: utf-8 -*-
using ProgressMeter
using FFTW, LinearAlgebra
using LinearMaps, IterativeSolvers
using Statistics
using Plots
using BenchmarkTools
using Test

# +
abstract type AbstractModel end
abstract type TimeSolver end
abstract type InitialData end

include("../src/data.jl")
include("../src/times.jl")
include("../src/mesh.jl")
include("../src/problem.jl")
include("../src/solvers/RK4.jl")
include("../src/initialdata/Init.jl")
include("../src/models/WhithamGreenNaghdi.jl")
# -

struct WhithamGreenNaghdiGPU <: AbstractModel

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

    function WhithamGreenNaghdiGPU(
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
            gtol
        )
    end
end


# +
"""
    mapto(WhithamGreenNaghdi, data)
`data` is of type `InitialData`, maybe constructed by `Init(...)`.

Performs a discrete Fourier transform with, possibly, dealiasing and Krasny filter.

See documentation of `WhithamGreenNaghdi` for more details.

"""
function mapto(m::WhithamGreenNaghdiGPU, data::InitialData)

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
function mapfro(m::WhithamGreenNaghdiGPU, datum::Array{Complex{Float64},2})
    m.fftη .= datum[:, 1]
    m.h .= 1 .+ m.ϵ * ifft(m.fftη)
    m.L .= m.Id .- 1 / 3 * m.FFT * Diagonal(1 ./ m.h) * m.M₀ * Diagonal(m.h .^ 3) * m.IFFTF₀

    real(ifft(datum[:, 1])), real(ifft(datum[:, 2])), real(ifft(m.L \ datum[:, 2]))
end


# +
function run_cpu( param )
    mesh = Mesh(param) # construct mesh of collocation points, Fourier modes, etc.
    η = exp.(-(mesh.x) .^ 2)
    v = 0 * η  # Initial data
    init = Init(mesh, η, v)
    model = WhithamGreenNaghdi(param) # iterate = true for GMRES iterate solver
    problem = Problem(model, init, param)
    U = copy(last(problem.data.U))
    @showprogress 1 for i in 1:param.ns
        step!(problem.solver, problem.model, U, problem.times.dt)
    end
    U
end

function run_gpu( param )
    mesh = Mesh(param) # construct mesh of collocation points, Fourier modes, etc.
    η = exp.(-(mesh.x) .^ 2)
    v = 0 * η  # Initial data
    init = Init(mesh, η, v)
    model = WhithamGreenNaghdiGPU(param) # iterate = true for GMRES iterate solver
    problem = Problem(model, init, param)
    U = copy(last(problem.data.U))
    @showprogress 1 for i in 1:param.ns
        step!(problem.solver, problem.model, U, problem.times.dt)
    end
    U
end
# -
function (m::WhithamGreenNaghdi)(U::Array{Complex{Float64},2})
    m.fftη .= U[:,1]
    m.h .= 1 .+ m.ϵ*ifft(m.fftη)
    m.fftv .= U[:,2]
    m.L .= m.Id - 1/3 * m.FFT * Diagonal( 1 ./m.h ) * m.M₀ * Diagonal( m.h.^3 ) * m.IFFTF₀
    m.fftu .= m.L \ m.fftv
    m.u .= ifft(m.fftu)
    m.hdu .= m.h .* ifft(m.F₀.*m.fftu)
    U[:,1] .= -m.∂ₓ.*fft(m.h .* m.u)
    U[:,2] .= -m.∂ₓ.*(m.fftη .+ m.ϵ * m.Π⅔.*fft( m.u.*ifft(m.fftv) .- 1/2 * m.u.^2 .- 1/2 * m.hdu.^2 ) )
    U[abs.(U).< m.ktol ].=0
end

# +
using LinearAlgebra.LAPACK: gesv!

function (m::WhithamGreenNaghdiGPU)(U::Array{Complex{Float64},2})
    m.fftη .= U[:, 1]
    m.h .= 1 .+ m.ϵ * ifft(m.fftη)
    m.fftv .= U[:, 2]

    m.L .= I - 1 / 3 .* m.FFT * Diagonal( 1 ./ m.h) * m.M₀ * Diagonal(m.h.^3) * m.IFFTF₀
    m.fftu .= m.fftv
    gesv!(m.L, m.fftu)
    m.u .= m.fftu
    ifft!(m.u)
    m.fftu .*= m.F₀
    ifft!(m.fftu)
    m.hdu .= m.h .* m.fftu
    
    m.h .*= m.u
    fft!(m.h)
    ifft!(m.fftv)
    m.u .= m.u .* m.fftv .- 1/2 .* m.u .^ 2 .- 1 / 2 .* m.hdu .^ 2
    fft!(m.u)
    U[:, 2] .= -m.∂ₓ .* ( view(U, :, 1) .+ m.ϵ * m.Π⅔ .* m.u )
    U[:, 1] .= -m.∂ₓ .* m.h
    
    U[abs.(U).<m.ktol] .= 0
end


# +
param = (
        μ = 1,
        ϵ = 1,
        N = 8,         # number of collocation points
        L = 10π,        # mesh of size 2L
        T = 1.0,        # final time
        dt = 0.0001,  # timestep
        ns = 1,      # stores data every ns time steps
    )

@test all( run_cpu(param) .== run_gpu(param))

# +
param = (
        μ = 1,
        ϵ = 1,
        N = 64,         # number of collocation points
        L = 10π,        # mesh of size 2L
        T = 1.0,        # final time
        dt = 0.0001,  # timestep
        ns = 10,      # stores data every ns time steps
    )

@time run_cpu(param)
@time run_gpu(param);
# -


