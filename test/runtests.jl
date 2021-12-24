ENV["GKSwstype"]="100"

using Test
using WaterWaves1D

include("../src/models/WaterWaves.jl")
include("../src/models/WWn.jl")

param = ( ϵ  = 1/2, μ = 1)
paramX= ( N  = 2^8, L  = 10)
paramT= ( T  = 5, dt = 0.1)

init     = Init(x->exp.(-x.^2),x->0 .*x )
model    = WaterWaves(merge(param,paramX);verbose=false)
problem1 = Problem( model, init, merge(paramT,paramX); solver = RK4(paramX) )

# @testset "LoadSave" begin
#
#     dump  = convert( ProblemSave, problem1 )
#     pload = convert( Problem, dump )
#
#     save(problem1, "testsave")
#     pload = loadpb("testsave")
#
#     @test pload.model.kwargs == problem1.model.kwargs
#     @test pload.solver.Uhat == problem1.solver.Uhat
#
# end

@testset "Parameters" begin

    @test param.ϵ  == 0.5
    @test param.μ  == 1
    @test paramX.N  == 256
    @test paramX.L  == 10
    @test paramT.T  == 5
    @test paramT.dt == 0.1

end

solve!(problem1)

@testset "Test problem with water waves model" begin
    @test !any(isnan,problem1.data.U[end][1])
    @test !any(isnan,problem1.data.U[end][2])
end

WW2  = WWn(merge(param,paramX);δ=0.01,n=2,dealias=1,verbose=false)
problem2 = Problem(WW2, init, merge(paramX,paramT) )

solve!( problem2 )

@testset "Test problem with Pseudo-spectral model" begin
    @test !any(isnan,problem2.data.U[end][1])
    @test !any(isnan,problem2.data.U[end][2])
end
