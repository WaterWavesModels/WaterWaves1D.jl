using Test

@testset "LoadSave" begin
 
    save_filename = "testsave"
    rm(joinpath(save_filename * ".h5"), force=true) # delete pre-existing file


    param  = ( c = 1.1, ϵ = 1, μ = 1, N  = 2^6, L = 4, T = 1, dt = 0.1)
    init = SolitaryWhitham(param)
    dump(save_filename, Mesh(param).x, init)

    model = Airy(param)
    problem = Problem( model, init, param )
    solve!(problem; verbose=false)

    dump(save_filename, problem.data)

    
    loaded_init = load_init(save_filename)
    @test loaded_init.label == "Whitham solitary wave"
    @test all( loaded_init.η(x) .≈ init.η(x) )
    @test all( loaded_init.v(x) .≈ init.v(x) )


    loaded_data = load_data(save_filename)
    @test loaded_data == problem.data

    rm(joinpath(save_filename * ".h5"), force=true) # delete pre-existing file

    dump(save_filename, problem)
    loaded_data = load_data(save_filename)
    @test loaded_data == problem.data

    new_problem = Problem( model, init, param )
    @test new_problem !== problem
    load_data!(save_filename, new_problem)
    @test new_problem == problem


end
