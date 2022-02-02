ENV["GKSwstype"]="100"

using Test
using WaterWaves1D

include("./test_loadsave.jl")
include("./testtypes.jl")
include("./testmodels.jl")
include("./testsolvers.jl")
include("./testinit.jl")


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

@testset "Interpolate" begin
    f(x)=exp.(-x.^2);
    mesh=Mesh(paramX);x=mesh.x;f1=f(x);
    mesh2,f2 = interpolate(mesh,f1;n=2)
    @test f2≈f(mesh2.x)
    @test f2[1:2:end]≈f1
    x3 = [1 2.5 3]
    f3 = interpolate(mesh,f1,x3)
    @test f3≈f(x3)
    x4 = [1 2 3]
    f4 = interpolate(mesh,f1,x4;fast=true)
    @test f4≈f(x4)
end
