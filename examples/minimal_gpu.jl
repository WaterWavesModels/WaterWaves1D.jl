using ProgressMeter
using FFTW, LinearAlgebra
using Test
using ShallowWaterModels

using CUDAdrv
using CuArrays
using Adapt

function run_cpu( param, steps )

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

    @showprogress 1 for i in 1:steps

        solver.Uhat .= U
        model( solver.Uhat )
        solver.Uhat .= U
        solver.dU .= solver.Uhat

        solver.Uhat .= U .+ dt/2 .* solver.Uhat
        model( solver.Uhat )
        solver.dU .+= 2 .* solver.Uhat

        solver.Uhat .= U .+ dt/2 .* solver.Uhat
        model( solver.Uhat )
        solver.dU .+= 2 .* solver.Uhat

        solver.Uhat .= U .+ dt .* solver.Uhat
        model( solver.Uhat )
        solver.dU .+= solver.Uhat

        U .+= dt/6 .* solver.dU

    end

    real(U)

end





function run_gpu( param, steps )

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

    d_h    = CuArray(model.h)
    d_u    = CuArray(model.u)
    d_hdu  = CuArray(model.hdu)
    d_fft  = CuArray(model.FFT)
    d_ifft = CuArray(model.IFFTF₀)
    d_mo   = CuArray(model.M₀)
    d_pi   = CuArray(model.Π⅔)
    d_id   = CuArray(model.Id)
    d_fo   = CuArray(model.F₀)
    d_L    = CuArray(model.L)
    d_dx   = CuArray(model.∂ₓ)

    function compute!(m::WhithamGreenNaghdi, U )
    
    	d_fftη = CuArray(U[:,1])
    	d_fftu = CuArray(U[:,2])
        d_fftv = copy(d_fftu)
    
        fw = CUFFT.plan_fft!(d_fftη)
        bw = CUFFT.plan_ifft!(d_fftv)
    
        d_h .= d_fftη
        bw * d_h
        d_h .*= m.ϵ
    	d_h .+= 1 
    
    	d_L .= Diagonal( d_h .* d_h .* d_h ) * d_ifft
    	d_L .= d_mo * d_L
    	d_L .= Diagonal( 1 ./ d_h ) * d_L
    	d_L .= d_fft * d_L
    	d_L .= d_id - 1/3 * d_L
    
        d_L, ipiv = CuArrays.CUSOLVER.getrf!(d_L)
        CuArrays.CUSOLVER.getrs!('N', d_L, ipiv, d_fftu)
    
        d_u .= d_fftu
        d_hdu .= d_u
        
        bw * d_u
    
    	d_hdu .*= d_fo
    	d_hdu .*= d_pi
        bw * d_hdu
    	d_hdu .*= d_h 
    	d_hdu .*= d_hdu
    	d_hdu .*= 1/2
    
        bw * d_fftv
        d_fftv .*= d_u
        d_fftv .-= 1/2 .* d_u .* d_u
        d_fftv .-= d_hdu
    
        fw * d_fftv
        d_fftv .*= m.ϵ
        d_fftv .+= d_fftη
        d_fftv .*= d_pi
        d_fftv .*= d_dx
    
        bw * d_fftη
        d_fftη .*= d_u
        fw * d_fftη
        d_fftη .*= m.ϵ
        d_fftη .+= d_fftu
        d_fftη .*= d_pi
        d_fftη .*= d_dx

        d_fftη .*= -1
        d_fftv .*= -1
    
    	U[:,1] .= Array(d_fftη)
    	U[:,2] .= Array(d_fftv)
    
    	U[abs.(U).< m.ktol ] .= 0.0
    
    end

    @showprogress 1 for i in 1:steps

        solver.Uhat .= U
        compute!(model, solver.Uhat )
        solver.Uhat .= U
        solver.dU .= solver.Uhat

        solver.Uhat .= U .+ dt/2 .* solver.Uhat
        compute!(model, solver.Uhat )
        solver.dU .+= 2 .* solver.Uhat

        solver.Uhat .= U .+ dt/2 .* solver.Uhat
        compute!(model, solver.Uhat )
        solver.dU .+= 2 .* solver.Uhat

        solver.Uhat .= U .+ dt .* solver.Uhat
        compute!(model, solver.Uhat )
        solver.dU .+= solver.Uhat

        U .+= dt/6 .* solver.dU

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

# trigger compilation
run_gpu(param, 1)
run_cpu(param, 1)

reset_timer!()

@show CUDAdrv.name(CuDevice(0))

@time U_gpu = run_gpu(param, 100)
@time U_cpu = run_cpu(param, 100)

@test U_cpu ≈ U_gpu
