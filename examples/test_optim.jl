using ProgressMeter
using FFTW, LinearAlgebra
using Plots
using JLD

abstract type AbstractModel end
abstract type TimeSolver    end
abstract type InitialData   end

mutable struct Data
    U :: Array{Array{Complex{Float64},2}}
    datasize :: Int
    datalength:: Int

    function Data( v )
        (datalength , datasize ) = size(v)
        U = Array{ComplexF64,2}[]
        push!(U, v)
        new(U, datasize, datalength)
    end
end

struct Mesh

    N    :: Int64
    xmin :: Float64
    xmax :: Float64
    dx   :: Float64
    x    :: Vector{Float64}
    kmin :: Float64
    kmax :: Float64
    dk   :: Float64
    k    :: Vector{Float64}

    function Mesh( xmin :: Float64, xmax :: Float64, N :: Int64)

        dx   = (xmax-xmin)/N
        x    = zeros(Float64, N)
        x   .= range(xmin, stop=xmax, length=N+1)[1:end-1] 
        dk   = 2π/(N*dx)
        kmin = -N/2*dk
        kmax = (N/2-1)*dk
        k    = zeros(Float64, N)
        k   .= dk .* vcat(0:N÷2-1, -N÷2:-1)

        new( N, xmin, xmax, dx, x, kmin, kmax, dk, k)

    end

    function Mesh(param :: NamedTuple)

        xmin = - Float64(param.L)
        xmax =   Float64(param.L)
        N    =   param.N
        
        Mesh( xmin, xmax, N)

    end
end

struct Times

    Nt   :: Int
    Nr   :: Int
    nr   :: Int
    tfin :: Float64
    dt   :: Float64
    t    :: Vector{Float64}
    tr   :: Vector{Float64}

    function Times( dt , tfin)
        t  = 0:dt:tfin
        Nt = length(t)
        Nr = Nt
        nr = 1
        tr = t
        new( Nt, Nr, nr, tfin, dt, t, tr)
    end

end


mutable struct Problem

    model   :: AbstractModel
    initial :: InitialData
    param   :: NamedTuple
    solver  :: TimeSolver
    times   :: Times
    mesh    :: Mesh
    data    :: Data

    function Problem( model   :: AbstractModel,
                      initial :: InitialData,
                      param   :: NamedTuple)

        times  = Times(param.dt, param.T)
        mesh   = Mesh(param)
        data   = Data(mapto(model,initial))
        solver = RK4(param,model)

        new(model,initial,param,solver,times,mesh,data)

    end

end

"""
    RK4(params)

Runge-Kutta fourth order solver.

"""
mutable struct RK4 <: TimeSolver

    Uhat :: Array{Complex{Float64},2}
    dU   :: Array{Complex{Float64},2}

    function RK4( param::NamedTuple, model::AbstractModel )

        n = param.N

        Uhat = zeros(Complex{Float64}, (n,model.datasize))
        dU   = zeros(Complex{Float64}, (n,model.datasize))

        new( Uhat, dU)

    end

end

function step!(s  :: RK4,
               f! :: AbstractModel,
               U  :: Array{Complex{Float64},2},
               dt :: Float64)

    s.Uhat .= U

    f!( s.Uhat )

    @simd for i in eachindex(U)
        @inbounds s.dU[i]   = s.Uhat[i]
        @inbounds s.Uhat[i] = U[i] + dt/2 * s.Uhat[i]
    end

    f!( s.Uhat )

    @simd for i in eachindex(U)
        @inbounds s.dU[i]   += 2 * s.Uhat[i]
        @inbounds s.Uhat[i]  = U[i] + dt/2 * s.Uhat[i]
    end

    f!( s.Uhat )

    @simd for i in eachindex(U)
        @inbounds s.dU[i]   += 2 * s.Uhat[i]
        @inbounds s.Uhat[i]  = U[i] + dt * s.Uhat[i]
    end

    f!( s.Uhat )

    @simd for i in eachindex(U)
        @inbounds s.dU[i] += s.Uhat[i]
        @inbounds U[i]    += dt/6 * s.dU[i]
    end

end

struct BellCurve <: InitialData

    h :: Vector{Float64}
    u :: Vector{Float64}

    function BellCurve(p :: NamedTuple,theta :: Real)

        mesh  = Mesh(p)
        h     = zeros(Float64, mesh.N)
        h    .= exp.(.-((abs.(mesh.x)).^theta).*log(2))
        u     = zeros(Float64, mesh.N)

        new( h, u )

    end

end

"""
    Matsuno(params)

"""
mutable struct Matsuno <: AbstractModel

    label    :: String
    datasize :: Int
    Γ        :: Array{Float64,1}
    Dx       :: Array{Complex{Float64},1}
    H        :: Array{Complex{Float64},1}
    Π⅔       :: BitArray{1}
    ϵ        :: Float64
    hnew     :: Vector{Complex{Float64}}
    unew     :: Vector{Complex{Float64}}
    I₀       :: Vector{Complex{Float64}}
    I₁       :: Vector{Complex{Float64}}
    I₂       :: Vector{Complex{Float64}}
    I₃       :: Vector{Complex{Float64}}

    Px       :: FFTW.FFTWPlan

    function Matsuno(param::NamedTuple)

        label    = "Matsuno"
        datasize = 2
        ϵ        = param.ϵ
        mesh     = Mesh(param)
        Γ        = abs.(mesh.k)
        Dx       =  1im * mesh.k        # Differentiation
        H        = -1im * sign.(mesh.k) # Hilbert transform
        Π⅔       = Γ .< mesh.kmax * 2/3 # Dealiasing low-pass filter

        hnew = zeros(Complex{Float64}, mesh.N)
        unew = zeros(Complex{Float64}, mesh.N)

        I₀ = zeros(Complex{Float64}, mesh.N)
        I₁ = zeros(Complex{Float64}, mesh.N)
        I₂ = zeros(Complex{Float64}, mesh.N)
        I₃ = zeros(Complex{Float64}, mesh.N)

        Px  = plan_fft(hnew; flags = FFTW.MEASURE)

        new(label, datasize, Γ, Dx, H, Π⅔, ϵ,
            hnew, unew, I₀, I₁, I₂, I₃, Px)
    end
end


function (m::Matsuno)(U::Array{Complex{Float64},2})


    @simd for i in eachindex(m.hnew)
        @inbounds m.hnew[i] = m.Γ[i] * U[i,1]
    end

    ldiv!(m.unew, m.Px, m.hnew )

    @simd for i in eachindex(m.hnew)
        @inbounds m.hnew[i] = m.Dx[i] * U[i,1]
    end

    ldiv!(m.I₁, m.Px, m.hnew)

    @simd for i in eachindex(m.unew)
        @inbounds m.unew[i] *= m.I₁[i]
    end

    mul!(m.I₁, m.Px, m.unew)

    @simd for i in eachindex(m.hnew)
        @inbounds m.I₁[i] = m.I₁[i] * m.ϵ[i] * m.Π⅔[i] - m.hnew[i]
    end

    ldiv!(m.hnew, m.Px, view(U,:,1))
    ldiv!(m.unew, m.Px, view(U,:,2))

    @simd for i in eachindex(m.hnew)
        @inbounds m.I₂[i] = m.hnew[i] * m.unew[i]
    end

    mul!(m.I₃, m.Px, m.I₂)

    @simd for i in eachindex(m.H)
        @inbounds U[i,1]  = m.H[i] * U[i,2]
        @inbounds m.I₀[i] = m.Γ[i] * U[i,2]
    end

    ldiv!(m.I₂, m.Px, m.I₀)

    @simd for i in eachindex(m.hnew)
        @inbounds m.I₂[i] *= m.hnew[i]
    end

    mul!(m.hnew, m.Px, m.I₂)

    @simd for i in eachindex(m.unew)
        @inbounds U[i,1] -= (m.I₃[i] * m.Dx[i] + m.hnew[i] * m.H[i]) * m.ϵ * m.Π⅔[i]
        @inbounds m.I₃[i]  = m.unew[i]^2
    end 

    mul!(m.unew, m.Px, m.I₃)

    @simd for i in eachindex(m.unew)
        @inbounds U[i,2]  =  m.I₁[i] - m.unew[i] * m.Dx[i] * m.ϵ/2 * m.Π⅔[i]
    end 

end

"""
    mapto(Matsuno, data)

"""
function mapto(m::Matsuno, data::InitialData)

    [m.Π⅔ .* fft(data.h) m.Π⅔ .* fft(data.u)]

end

"""
    mapfro(Matsuno, data)

"""
function mapfro(m::Matsuno, datum::Array{Complex{Float64},2})

    real(ifft(view(datum,:,1))),real(ifft(view(datum,:,2)))

end

function create_animation( p::Problem )

    prog = Progress(p.times.Nr,1)

    anim = @animate for (l,U) in enumerate(p.data.U)

        pl = plot(layout=(2,1))

		(hr,ur) = mapfro(p.model,U)

        plot!(pl[1,1], p.mesh.x, hr;
              ylims=(-0.6,1),
              title="physical space",
              label=p.model.label)

        plot!(pl[2,1], fftshift(p.mesh.k),
              log10.(1e-18.+abs.(fftshift(fft(hr))));
              title="frequency",
              label=p.model.label)

        next!(prog)

    end when mod(l, 100) == 0

    gif(anim, "anim.gif", fps=15); nothing

end

function solve!(problem :: Problem)

    @show problem.param

    U  = similar(problem.solver.Uhat)
    U .= problem.data.U[end]
    step!(problem.solver, problem.model, U, 0.0)

    dt = problem.times.dt

    nr = problem.times.nr

    J = nr:nr:problem.times.Nt-1
    L = 1:nr

    @showprogress 1 for j in J
        for l in L
            step!(problem.solver, problem.model, U, dt)
        end
        push!(problem.data.U,copy(U))
    end

    print("\n")

end


function main()

    param = ( ϵ  = 1/2,
              N  = 2^12,
              L  = 10.,
              T  = 5.,
              dt = 0.001 )
    
    init    = BellCurve(param,2.5)
    model   = Matsuno(param)
    problem = Problem(model, init, param);
    
    @time solve!( problem )

    Uref =  problem.data.U[end]

    #save("reference.jld", "Uref", Uref)

    Uref = load("reference.jld", "Uref")

    println(norm(Uref .- problem.data.U[end]))

    create_animation( problem )

end

main()
