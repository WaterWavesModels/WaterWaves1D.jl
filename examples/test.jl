using ProgressMeter
using FFTW, LinearAlgebra
using Plots
gr()

abstract type AbstractModel end
abstract type TimeSolver end
abstract type InitialData end

mutable struct Data
    U :: Array{Array{Complex{Float64},2}}
    datasize :: Int
    datalength:: Int

    function Data( v )
        (datalength , datasize ) = size(v)
        U = [v]
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

    function Times( dt, tfin)
        t = range(0, stop=tfin, step = dt)
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
    s.dU .= s.Uhat

    s.Uhat .= U .+ dt/2 .* s.Uhat
    f!( s.Uhat )
    s.dU .+= 2 .* s.Uhat

    s.Uhat .= U .+ dt/2 .* s.Uhat
    f!( s.Uhat )
    s.dU .+= 2 .* s.Uhat

    s.Uhat .= U .+ dt .* s.Uhat
    f!( s.Uhat )
    s.dU .+= s.Uhat

    U .+= dt/6 .* s.dU

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


    for i in eachindex(m.hnew)
        m.hnew[i] = m.Γ[i] * U[i,1]
    end

    ldiv!(m.unew, m.Px, m.hnew )

    for i in eachindex(m.hnew)
        m.hnew[i] = m.Dx[i] * U[i,1]
    end

    ldiv!(m.I₁, m.Px, m.hnew)

    m.unew  .*= m.I₁

    mul!(m.I₁, m.Px, m.unew)

    m.I₁  .*= m.ϵ .* m.Π⅔
    m.I₁  .-= m.hnew

    ldiv!(m.hnew, m.Px, view(U,:,1))
    ldiv!(m.unew, m.Px, view(U,:,2))

    m.I₂    .= m.hnew .* m.unew

    mul!(m.I₃, m.Px, m.I₂)

    m.I₃    .*= m.Dx

    for i in eachindex(m.H)
        U[i,1]  = m.H[i] * U[i,2]
        m.I₀[i] = m.Γ[i] * U[i,2]
    end

    ldiv!(m.I₂, m.Px, m.I₀)

    m.I₂    .*= m.hnew

    mul!(m.hnew, m.Px, m.I₂)

    m.hnew  .*= m.H
    m.I₃    .+= m.hnew
    m.I₃    .*= m.ϵ .* m.Π⅔
    
    for i in eachindex(m.I₃)
        U[i,1] -= m.I₃[i]
    end 

    m.I₃    .=  m.unew.^2

    mul!(m.unew, m.Px, m.I₃)

    m.unew  .*= m.Dx
    m.unew  .*= m.ϵ/2 .* m.Π⅔
    m.I₁    .-= m.unew

    for i in eachindex(m.I₁)
        U[i,2] =  m.I₁[i]
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


function solve!(problem :: Problem)

    @show problem.param

    U  = similar(problem.solver.Uhat)
    U .= problem.data.U[end]
    step!(problem.solver, problem.model, U, 0.0)

    dt = problem.times.dt

    prog = Progress(problem.times.Nt,1)

    J = range(problem.times.nr ,stop = problem.times.Nt-1, step = problem.times.nr)
    L = 1:problem.times.nr

    for j in J
        for l in L

            step!(problem.solver, problem.model, U, dt)

            next!(prog)

        end

        push!(problem.data.U,copy(U))


    end

    print("\n")

end


function main()

    param = ( ϵ  = 1/2,
              N  = 2^12,
              L  = 10,
              T  = 5,
              dt = 0.001 )
    
    init    = BellCurve(param,2.5)
    model   = Matsuno(param)
    problem = Problem(model, init, param);
    
    solve!( problem )

   # (hr,ur) = mapfro(problem.model,problem.data.U[end])

   # p = plot(layout=(2,1))
   # 
   # plot!(p[1,1], problem.mesh.x, hr;
   #     	  title="physical space",
   #               label=problem.model.label)

   # plot!(p[2,1], fftshift(problem.mesh.k),
   #               log10.(1e-18.+abs.(fftshift(fft(hr))));
   #     	  title="frequency",
   # 	          label=problem.model.label)

end

@time main()
