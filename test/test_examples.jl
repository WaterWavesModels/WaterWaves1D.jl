using Test
import LinearAlgebra.norm
@testset "Example: ManyModelsforWaterWaves" begin
    include("../examples/ManyModelsforWaterWaves.jl")

    Diff,=Integrate(1;μ=.1,ϵ=.1,p=2,N=2^8,L=10,T=1,dt=0.01,
            dealias=0,iterate=true,precond=true,method=1,maxiter=10,
            name="test_1")
    
    @test 1e-3< Diff[1] <2e-3   # SGN
    @test 1e-4< Diff[2] <4e-4   # WGN
    @test 1e-5< Diff[3] <4e-5   # IK2

    Diff,=Integrate(2;μ=.1,ϵ=.1,p=2,N=2^8,L=10,T=1,dt=0.01,
            dealias=0,iterate=true,precond=true,method=2,maxiter=10,
            name="test_2")
            
    @test 2e-3< Diff[1] <4e-3   # SGN
    @test 4e-4< Diff[2] <8e-4   # WGN
    @test 1e-4< Diff[3] <2e-4   # IK2
        

end
@testset "Example: ManyModelsforWaterWaves" begin
    include("../examples/SolitaryWaveWhitham.jl")

    a=PlotSolitaryWaveKdV(3)
    @test a ≈ 4.5488501854151764e-11

    b=PlotSolitaryWaveWhitham(1.17)
    @test b ≈  0.06741522566084218

end