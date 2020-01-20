ENV["GKSwstype"]="100"

using Test
using ShallowWaterModels
using JLD

param = ( ϵ  = 1/2,
          N  = 2^10,
          L  = 10,
          T  = 5,
          dt = 0.01,
          theta = 1)

init     = BellCurve(param)
solver   = RK4(param)
cheng    = CGBSW(param)
problem1 = Problem( cheng, init, param, solver )

@testset "LoadSave" begin

    dump  = convert(ProblemSave, problem1 )
    p2    = convert(Problem, dump )

    save(problem1, "problem1")

    pload = load("problem1")

    @test true

end

@testset "Parameters" begin

    @test param.ϵ  == 0.5
    @test param.N  == 1024
    @test param.L  == 10
    @test param.T  == 5
    @test param.dt == 0.01
    @test param.theta == 1

end

solve!(problem1)

@testset "Test problem with Cheng et al. model" begin
    @test !any(isnan,problem1.data.U[end][1])
    @test !any(isnan,problem1.data.U[end][2])
end

matsuno  = Matsuno(param)
problem2 = Problem(matsuno, init, param, solver )

solve!( problem2 )

@testset "Test problem with Matsuno model" begin
    @test !any(isnan,problem2.data.U[end][1])
    @test !any(isnan,problem2.data.U[end][2])
end
