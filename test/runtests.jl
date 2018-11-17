using Test
using DeepWaterModels
using FFTW

param = Parameters( ϵ  = 1/2,
                    N  = 2^10,
                    L  = 10,
                    T  = 5,
                    dt = 0.01)


@testset "Parameters" begin

    @test param.ϵ  == 0.5
    @test param.N  == 1024
    @test param.L  == 10
    @test param.T  == 5
    @test param.dt == 0.01

end

init     = BellCurve(param,1)
solver   = RK4(param)
cheng    = CGBSW(param)
problem1 = Problem( cheng, init, param, solver )

solve!(problem1)

@testset "Test problem with Cheng et al. model" begin
    @test !any(isnan,problem1.data[end][1])
    @test !any(isnan,problem1.data[end][2])
end

matsuno  = Matsuno(param)
problem2 = Problem(matsuno, init, param, solver )

solve!( problem2 )

@testset "Test problem with Matsuno model" begin
    @test !any(isnan,problem2.data[end][1])
    @test !any(isnan,problem2.data[end][2])
end
