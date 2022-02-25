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
paramT= ( T  = 0.1, dt = 0.01)

@testset "Parameters" begin

    @test param.ϵ  == 0.5
    @test param.μ  == 1
    @test paramX.N  == 256
    @test paramX.L  == 10
    @test paramT.T  == 0.1
    @test paramT.dt == 0.01

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

models=[];
push!(models,Airy(param;mesh=Mesh(paramX)));
push!(models,SaintVenant(param;mesh=Mesh(paramX)));
push!(models,Boussinesq(param;mesh=Mesh(paramX)));
push!(models,WhithamBoussinesq(param;mesh=Mesh(paramX)));
push!(models,SerreGreenNaghdi(param;mesh=Mesh(paramX)));
push!(models,WhithamGreenNaghdi(param;mesh=Mesh(paramX)));
push!(models,NonHydrostatic(param;mesh=Mesh(paramX)));
push!(models,SquareRootDepth(param;mesh=Mesh(paramX)));
push!(models,IsobeKakinuma(param;mesh=Mesh(paramX)));
push!(models,WWn(param;mesh=Mesh(paramX)));

for model in models
    @testset "Preserved quantities : $(model.label)" begin
    pb=Problem(model,Init(x->exp.(-x.^2),x->(x.+1).*exp.(-(x.+1).^2)),paramT)
    solve!(pb;verbose=false)
    @test isapprox(mass(pb),mass(pb;t=0),rtol=1e-10)
    @test isapprox(momentum(pb),momentum(pb;t=0),rtol=1e-10)
    @test isapprox(energy(pb),energy(pb;t=0),rtol=1e-8)
    @test abs(massdiff(pb))<1e-10
    @test abs(momentumdiff(pb))<1e-10
    @test abs(energydiff(pb;rel=true))<1e-8
    end
end
