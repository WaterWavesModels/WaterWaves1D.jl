using TimerOutputs
using ProgressMeter
using FFTW, LinearAlgebra
using Statistics
using JLD2, FileIO
using Printf
using Test
using ShallowWaterModels

function run_vincent( param )

    @timeit "Setup" begin

        mesh = Mesh(param) # construct mesh of collocation points, Fourier modes, etc.

        η = exp.(-(mesh.x) .^ 2)

        v = zero(η)  # Initial data

        init = Init(mesh, η, v)

        times = Times(param.dt, param.T)

	    model = WhithamGreenNaghdi(param; iterate=false)

        solver = RK4(param, model) 

        data  = Data(mapto(model, init))

        U = copy(last(data.U))

        dt = times.dt
    end

    @showprogress 1 for i in 1:param.ns

        @timeit "RK4" solver.Uhat .= U
        @timeit "Model" model( solver.Uhat )
        @timeit "RK4" solver.Uhat .= U
        @timeit "RK4" solver.dU .= solver.Uhat

        @timeit "RK4" solver.Uhat .= U .+ dt/2 .* solver.Uhat
        @timeit "Model" model( solver.Uhat )
        @timeit "RK4" solver.dU .+= 2 .* solver.Uhat

        @timeit "RK4" solver.Uhat .= U .+ dt/2 .* solver.Uhat
        @timeit "Model" model( solver.Uhat )
        @timeit "RK4" solver.dU .+= 2 .* solver.Uhat

        @timeit "RK4" solver.Uhat .= U .+ dt .* solver.Uhat
        @timeit "Model" model( solver.Uhat )
        @timeit "RK4" solver.dU .+= solver.Uhat

        @timeit "RK4" U .+= dt/6 .* solver.dU

        @timeit "Diagnostic" push!(data.U, copy(U))

    end

    real(U)

end



function compute!(m::WhithamGreenNaghdi, U )

	m.h .= U[:,1]
    ifft!(m.h)
    m.h .*= m.ϵ
	m.h .+= 1 

    
	@timeit "LU" m.L .= m.Id - 1/3 * m.FFT * Diagonal( 1 ./m.h ) * m.M₀ * Diagonal( m.h.^3 ) * m.IFFTF₀

	m.fftu .= m.L \ view(U,:,2)

	m.u .= ifft(m.fftu)

	m.hdu .= m.fftu
	m.hdu .*= m.F₀
	m.hdu .*= m.Π⅔
    ifft!(m.hdu)
	m.hdu .*= m.h 
	m.hdu .*= m.hdu
	m.hdu .*= 1/2

    m.fftv = U[:,2]
    ifft!(m.fftv)
    m.fftv .*= m.u
    m.fftv .-= 1/2 .* m.u.^2
    m.fftv .-= m.hdu
    fft!(m.fftv)
    m.fftv .*= m.ϵ
    m.fftv .+= view(U,:,1)
    m.fftv .*= m.Π⅔ 
    m.fftv .*= m.∂ₓ
	U[:,2] .= - m.fftv 


    m.fftv .= U[:,1]
    ifft!(m.fftv)
    m.fftv .*= m.u
    fft!(m.fftv)
    m.fftv .*= m.ϵ
    m.fftv .+= m.fftu
    m.fftv .*= m.Π⅔ 
    m.fftv .*= m.∂ₓ
   	U[:,1] .= - m.fftv


	U[abs.(U).< m.ktol ] .= 0.0

end


function run_pierre( param )

    @timeit "Setup" begin

        mesh = Mesh(param) # construct mesh of collocation points, Fourier modes, etc.

        η = exp.(-(mesh.x) .^ 2)

        v = zero(η)  # Initial data

        init = Init(mesh, η, v)

        times = Times(param.dt, param.T)

	    model = WhithamGreenNaghdi(param; iterate=false)

        solver = RK4(param, model) 

        data  = Data(mapto(model, init))

        U = copy(last(data.U))

        dt = times.dt
    end

    @showprogress 1 for i in 1:param.ns

        @timeit "RK4" solver.Uhat .= U
        @timeit "Model" compute!(model, solver.Uhat )
        @timeit "RK4" solver.Uhat .= U
        @timeit "RK4" solver.dU .= solver.Uhat

        @timeit "RK4" solver.Uhat .= U .+ dt/2 .* solver.Uhat
        @timeit "Model" compute!(model, solver.Uhat )
        @timeit "RK4" solver.dU .+= 2 .* solver.Uhat

        @timeit "RK4" solver.Uhat .= U .+ dt/2 .* solver.Uhat
        @timeit "Model" compute!(model, solver.Uhat )
        @timeit "RK4" solver.dU .+= 2 .* solver.Uhat

        @timeit "RK4" solver.Uhat .= U .+ dt .* solver.Uhat
        @timeit "Model" compute!(model, solver.Uhat )
        @timeit "RK4" solver.dU .+= solver.Uhat

        @timeit "RK4" U .+= dt/6 .* solver.dU

        @timeit "Diagnostic" push!(data.U, copy(U))

    end

    real(U)

end

param = (
        μ = 1,
        ϵ = 1,
        N = 1024,       # number of collocation points
        L = 10π,        # mesh of size 2L
        T = 1.0,        # final time
        dt = 0.0001,    # timestep
        ns = 10,         # stores data every ns time steps
      )

reset_timer!()

@time U_vincent = run_vincent(param)

print_timer()

reset_timer!()

@time U_pierre = run_pierre(param)

print_timer()

println()

@test U_vincent ≈ U_pierre
