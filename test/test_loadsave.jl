using Test
using WaterWaves1D
import WaterWaves1D: dump
import WaterWaves1D: load_data
import WaterWaves1D: load_mesh
import WaterWaves1D: load_times
import WaterWaves1D: load_timesolver
import WaterWaves1D: load_param
import WaterWaves1D: load_model
import WaterWaves1D: load_init

@testset "LoadSave" begin
 
    init  = Init(x->exp.(-1/8*x.^4),x-> 0. * x )

    para  = ( ϵ  = 0.01, μ = 0.01)  # physical parameters
    paraX = ( N  = 2^6, L  = 4)   # mesh with 64 collocation points on [-4,4]
    paraT = ( T  = 1e-1, dt = 1e-2) # timegrid with 10 instants: t=[0.0:1.0:10.0]/100
    param = merge(para,paraX)  # used to construct models
    parap = merge(paraX,paraT) # used to construct problems

    model = SerreGreenNaghdi(param;
            dealias = 0,
            ktol	= 1e-10,
            iterate = true,
            gtol	= 1e-14,
            precond = true,
            restart	= 15,
            maxiter	= 25,
			label	= "Green-Naghdi with GMRES")

    precision = para.μ^2

    problem = Problem( model, init, parap )
    solve!(problem; verbose=false)

    save_filename = "testsave"
    rm(joinpath(save_filename * ".h5"), force=true)

    dump(save_filename, model)
    @test load_model(save_filename) == "Green-Naghdi with GMRES"

    dump(save_filename, init)
    @test load_init(save_filename) == "user-defined"

    dump(save_filename, param)
    saved_param = load_param( save_filename )
    @test saved_param ≈ param

    dump(save_filename, problem.solver)
    solver = load_timesolver(save_filename)
    @test solver == problem.solver

    dump(save_filename, problem.times)
    times = load_times(save_filename)
    @test times == problem.times

    dump(save_filename, problem.mesh)
    mesh = load_mesh(save_filename)
    @test mesh == problem.mesh

    dump(save_filename, problem.data)
    data = load_data(save_filename)
    @test length(data) == length(problem.data)

 
end
