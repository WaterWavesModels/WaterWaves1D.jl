using Test
import LinearAlgebra.norm
@testset "Example: ManyModelsforWaterWaves" begin
    include("../examples/ManyModelsforWaterWaves.jl")

    Diff,=Integrate(1;μ=.1,ϵ=.1,p=2,N=2^8,L=10,T=1,dt=0.01,
            dealias=0,iterate=true,precond=true,method=1,maxiter=10,
            name=nothing)
    
    @test 1e-3< Diff[1] <2e-3   # SGN
    @test 1e-4< Diff[2] <4e-4   # WGN
    @test 1e-5< Diff[3] <4e-5   # IK2

    Diff,=Integrate(2;μ=.1,ϵ=.1,p=2,N=2^8,L=10,T=1,dt=0.01,
            dealias=0,iterate=true,precond=true,method=2,maxiter=10,
            name=nothing)
            
    @test 2e-3< Diff[1] <4e-3   # SGN
    @test 4e-4< Diff[2] <8e-4   # WGN
    @test 1e-4< Diff[3] <2e-4   # IK2
        

end

@testset "Example: Solitary Wave Whitham and KdV" begin
    include("../examples/SolitaryWaveWhitham.jl")

    a=PlotSolitaryWaveKdV(3)
    @test a ≈ 4.5488501854151764e-11

    b=PlotSolitaryWaveWhitham(1.17)
    @test b ≈  0.06741522566084218

end

@testset "Example: Solitary Wave Whitham-Boussinesq" begin
        include("../examples/SolitaryWaveWhithamBoussinesq.jl")

        b=IntegrateSolitaryWaveWhithamBoussinesq(N=2^6,dt=0.01/1.05)
        @test b ≈  9.904591036224986e-11
    
        a=PlotSolitaryWaveWhithamBoussinesq(;c=1.1,α=1,L=20,N=2^9,μ=0.1,ϵ=0.1)
        @test a ≈ 0.031638234910189356
    
    end