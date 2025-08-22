@testset "Interpolate" begin
    paramX= ( N  = 2^8, L  = 10)
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

@testset "Solution" begin

    param = ( ϵ  = 1/2, μ = 1, N  = 2^8, L  = 10 , T  = 0.1, dt = 0.01)
    init = Init(x->sin.(-x.^2),x->cos.(x)) 
    model = Airy(param;mesh=Mesh(param))
    pb = Problem(model,init,param)
    solve!(pb)

    η₀,v₀,x₀,t₀=solution(pb,T=0)
    @test η₀[1]≈sin(-x₀[1]^2)
    @test v₀[1]≈cos(x₀[1])
    
    η,v,x,t=solution(pb)
    @test x==Mesh(param).x
    @test t==param.T
    η₁,v₁,x₁,t₁=solution(pb;T=param.T,x=[x[1] 2])
    @test t₁==param.T
    @test x₁==[x[1] 2]
    @test η[1] ≈ η₁[1]
    @test v[1] ≈ v₁[1]

    η₂,v₂,x₂=solution(pb;T=param.T,interpolation=true)
    @test x₂[1:2^3:end] ≈ x
    @test η₂[1:2^3:end] ≈ η
    @test v₂[1:2^3:end] ≈ v

    η₃,v₃,x₃=solution(pb;T=param.T,interpolation=4)
    @test x₃[1:4:end] ≈ x
    @test η₃[1:4:end] ≈ η
    @test v₃[1:4:end] ≈ v

    raw_η,raw_v,t=solution(pb;raw=true)	
    @test raw_η==pb.data.U[end][1]
    @test raw_v==pb.data.U[end][2]
    @test t==pb.times.ts[end]
end

using Plots
@testset "Plot recipes" begin
    param = ( ϵ  = 1/2, μ = 1, N  = 2^8, L  = 10 , T  = 0.1, dt = 0.01)
    init = Init(x->sin.(-x.^2),x->cos.(x)) 
    model = Airy(param;mesh=Mesh(param))
    model2 = WWn(param;mesh=Mesh(param))
    model3 = WWn(param;mesh=Mesh(param),n=1)

    pb = Problem(model,init,param)
    pb2 = Problem(model2,init,param)
    pb3 = Problem(model3,init,param)

    solve!([pb,pb2,pb3])

    plt=plot(pb)
    @test plt.n==1

    plt=plot(pb,var=:surface)
    @test plt.n==1

    plt=plot(pb,var=[:surface,:velocity,:fourier,:Fourier,:fourier_surface , :Fourier_surface, :fourier_velocity , :Fourier_velocity];interpolation=true)
    @test plt.n==8

    plt=plot([pb,pb2],var=[:surface,:velocity,:fourier,:Fourier,:fourier_surface , :Fourier_surface, :fourier_velocity , :Fourier_velocity];compression=true)
    @test plt.n==16

    plt=plot([pb,pb2],var=[:difference , :difference_surface, :difference_velocity, :difference_fourier , :difference_Fourier])
    @test plt.n==5

    plt=plot([pb,pb2,pb3],var=[:differences , :differences_surface, :differences_velocity, :differences_fourier , :differences_Fourier])
    @test plt.n==15

    plt=plot((pb,pb2),var=[:surface,:velocity,:fourier])
    @test plt.n==3

    plt=plot([(pb,pb2),(pb,pb3)],var=[:surface_difference,:velocity_difference,:fourier_difference])
    @test plt.n==6

end

### Preserved quantities
param = ( ϵ  = 1/2, μ = 1)
paramX= ( N  = 2^8, L  = 10)
paramT= ( T  = 0.1, dt = 0.01)

models=AbstractModel[];
#push!(models,WaterWaves(param;mesh=Mesh(paramX))); # will error if uncommented (because water waves generates non-equidistant meshes)
#push!(models,Matsuno(param;mesh=Mesh(paramX))); # will fail if uncommented (because Matsuno model does not preserve momentum and energy)

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
    @test isapprox(mass(pb),mass(pb;T=0),rtol=1e-10)
    @test isapprox(momentum(pb),momentum(pb;T=0),rtol=1e-10)
    @test isapprox(energy(pb),energy(pb;T=0),rtol=1e-8)
    @test abs(massdiff(pb))<1e-10
    @test abs(momentumdiff(pb))<1e-10
    @test abs(energydiff(pb;rel=true))<1e-8
    end
end
