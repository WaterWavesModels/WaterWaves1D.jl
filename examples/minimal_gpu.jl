using ProgressMeter
using FFTW, LinearAlgebra
using Test
using ShallowWaterModels

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

	problem = Problem(model, init, param)

	solve!( problem )

    real(problem.U)

end

function run_gpu( param, steps )

    mesh = Mesh(param) # construct mesh of collocation points, Fourier modes, etc.
    η = exp.(-(mesh.x) .^ 2)
    v = zero(η)  # Initial data
    init = Init(mesh, η, v)
    times = Times(param.dt, param.T)
	model = WhithamGreenNaghdiGPU(param; iterate=false)
    solver = RK4(param, model) 
    data  = Data(mapto(model, init))
    U = copy(last(data.U))
    dt = times.dt

	problem = Problem(model, init, param)

	solve!( problem )

    real(problem.U)

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

@show CUDAdrv.name(CuDevice(0))

@time U_gpu = run_gpu(param, 100)
@time U_cpu = run_cpu(param, 100)

@test U_cpu ≈ U_gpu
