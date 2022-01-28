ENV["GKSwstype"]="100"

using Test
using WaterWaves1D

include("./testtypes.jl")
include("./testmodels.jl")
include("./testsolvers.jl")


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

param = ( ϵ  = 1/2, μ = 1)
paramX= ( N  = 2^8, L  = 10)
paramT= ( T  = 5, dt = 0.1)


@testset "Parameters" begin

    @test param.ϵ  == 0.5
    @test param.μ  == 1
    @test paramX.N  == 256
    @test paramX.L  == 10
    @test paramT.T  == 5
    @test paramT.dt == 0.1

end
