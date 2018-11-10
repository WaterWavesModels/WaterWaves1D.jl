using Test
using DeepWaterModels
using FFTW

param = Parameters( ϵ  = 1/2, 
                    N  = 2^12,
                    L  = 10,
                    T  = 5,
                    dt = 0.001)


@testset "Parameters" begin

    @test param.ϵ  == 0.5
    @test param.N  == 4096
    @test param.L  == 10
    @test param.T  == 5
    @test param.dt == 0.001

end

bump     = Bump(param)
solver   = RK4(param)
cheng    = Cheng(param)
problem1 = Problem( cheng, bump, param, solver )

solve!(problem1)

@testset "Test problem with Cheng model" begin
    @test !any(isnan,problem1.data[end][1])
    @test !any(isnan,problem1.data[end][2])
end

matsuno  = Matsuno(param)
problem2 = Problem(matsuno, bump, param, solver )

solve!( problem2 )

@testset "Test problem with Matsuno model" begin
    @test !any(isnan,problem2.data[end][1])
    @test !any(isnan,problem2.data[end][2])
end
