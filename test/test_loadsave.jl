using Test

@testset "Tests on LoadSave" begin
 
    save_filename = "testsave"
    rm(joinpath(save_filename * ".h5"), force=true) # delete pre-existing file


    param  = ( h₀ = 1, h₁ = 1.5, h₂ = 2, ϵ = 1, μ = 1, N  = 2^6)
    (η,u,v,mesh,par) = CnoidalWaveSerreGreenNaghdi(param)
    param=merge(par,param,( L=par.λ, T = 1, dt = 0.1))
    init = Init(mesh,η,v;label = "SGN cnoidal wave")
    x = mesh.x
    dump(save_filename, x, init)

    model = Airy(param)
    problem = Problem( model, init, param )
    solve!(problem; verbose=false)

    dump(save_filename, problem.data)

    
    loaded_init = load_init(save_filename)
    @test loaded_init.label == "SGN cnoidal wave"
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
