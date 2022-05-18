using Test

@testset "LoadSave" begin
 
    η = x -> exp.(-1/8*x.^4);
    v = zero
    init  = Init( η , v )

    param  = ( μ = 0.01, N  = 2^6, L  = 4, T  = 1e-1, dt = 1e-2)

    model = Airy(param)

    problem = Problem( model, init, param )
    solve!(problem; verbose=false)

    save_filename = "testsave"
    rm(joinpath(save_filename * ".h5"), force=true) # delete pre-existing file

    x = Mesh(param).x
    dump(save_filename, Mesh(param).x, init)
    loaded_init = load_init(save_filename)
    @test loaded_init.label == "user-defined"
    @test all(abs.(loaded_init.η(x)-η(x)) .< 4*eps())
    @test all(abs.(loaded_init.v(x)-v(x)) .< 4*eps())


    dump(save_filename, problem.data)
    data = load_data(save_filename)
    @test data.datalength == param.N
    @test data.datasize == 2
    @test data.U == problem.data.U

    rm(joinpath(save_filename * ".h5"), force=true) # delete pre-existing file

    dump(save_filename, problem)
    data = load_data(save_filename)
    @test data.datalength == param.N
    @test data.datasize == 2
    @test data.U == problem.data.U

    problem_loaded = Problem( model, init, param )
    @test problem_loaded !== problem
    load_data!(problem_loaded, save_filename)
    @test problem_loaded == problem




end
