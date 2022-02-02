using Test



#--- tests on SGN cnoidal wave
@testset "SGN cnoidal wave" begin
    # Build in a first way a cnoidal wave
    param = (ϵ  = 0.1, μ = 0.1, h₀ = 0.5, h₁ = 0.8, h₂ = 1.5, N=2^9)
    (η,u,v,mesh,para)=CnoidalWaveSerreGreenNaghdi(param;P=1)
    # Build the cnoidal wave with two periods
    (η2,u2,v2,mesh2,para2)=CnoidalWaveSerreGreenNaghdi(param;P=2)
    # Build in a second way the cnoidal wave
    init = CnoidalSGN(param)
    # Check they correspond
    @test init.η(mesh.x)≈η
    @test init.v(mesh.x)≈v
    @test init.η(mesh2.x)≈η2
    @test init.v(mesh2.x)≈v2
    @test η2[end÷2+1:3*end÷4]≈η[end÷2+1:2:end]
    @test v2[end÷2+1:3*end÷4]≈v[end÷2+1:2:end]

    # Solve the corresponding initial-value problem
    param=merge(param,(L=para.λ,T=1e-2,dt=1e-3))
    model = SerreGreenNaghdi(param)
    problem = Problem(model, init, param)
    solve!(problem;verbose=false)
    # Check the solution is translated
    (η,v,x,t)=solution(problem)
    @test η≈init.η(x.-para.c*param.T)
    @test v≈init.v(x.-para.c*param.T)

end

#--- tests on SGN solitary wave
@testset "SGN solitary wave" begin
    # Build in a first way a solitary wave
    param = (ϵ  = 0.1, μ = 0.1, c=1.1, L = 10, N=2^8)
    (η,u,v,mesh)=SolitaryWaveSerreGreenNaghdi(param)
    # Build the translated solitary wave on a different mesh
    param2 = (ϵ  = 0.1, μ = 0.1, c=1.1, L = 20, N=2^9);x₀=-2param.L/param.N;
    (η2,u2,v2,mesh2)=SolitaryWaveSerreGreenNaghdi(param2;x₀=x₀)
    # Build in a second way the cnoidal wave
    init = SolitarySGN(param2)
    init2 = SolitarySGN(param2.c,ϵ=param2.ϵ,μ=param2.μ,N=2^10)

    # Check they correspond
    @test init.η(mesh.x)≈η
    @test init.v(mesh.x)≈v
    @test init.η(mesh.x)≈init2.η(mesh.x)
    @test init.v(mesh.x)≈init2.v(mesh.x)
    @test init.η(mesh2.x .-x₀)≈η2
    @test init.v(mesh2.x .-x₀)≈v2
    @test η2[end÷4:3*end÷4-1]≈η[1:end]
    @test u2[end÷4:3*end÷4-1]≈u[1:end]
    @test v2[end÷4:3*end÷4-1]≈v[1:end]

    # Solve the corresponding initial-value problem
    param=merge(param2,(T=1e-2,dt=1e-3))
    model = SerreGreenNaghdi(param)
    problem = Problem(model, init, param)
    solve!(problem;verbose=false)
    # Check the solution is translated
    (η,v,x,t)=solution(problem)
    @test η≈init.η(x.-param.c*param.T)
    @test v≈init.v(x.-param.c*param.T)

end

#--- tests on Whitham-Boussinesq solitary wave
@testset "Whitham-Boussinesq solitary wave" begin
    # Build a solitary wave with different parameters
    param = (ϵ  = 0.1, μ = 0.1, c=1.1, L = 10, N=2^8)
    waves = []
    push!(waves,SolitaryWaveWhithamBoussinesq(param))
    push!(waves,SolitaryWaveWhithamBoussinesq(param;iterative=true))
    push!(waves,SolitaryWaveWhithamBoussinesq(param;
            iterative=true,max_iter = 15, tol = 1e-9, gtol = 1e-9, ktol = 1e-9))
    push!(waves,SolitaryWaveWhithamBoussinesq(param;dealias=1))
    push!(waves,SolitaryWaveWhithamBoussinesq(param;q=0.9,β=1/2,max_iter=30))
    push!(waves,SolitaryWaveWhithamBoussinesq(param;β=1/2))

    # Check they are consistent
    (η,v,mesh)=waves[1]
    for wave in waves
        (η2,v2,mesh2)=wave
        @test η2≈η
        @test v2≈v
    end

    # Build the translated solitary wave on a different mesh
    param2 = (ϵ  = 0.1, μ = 0.1, c=1.1, L = 20, N=2^9);x₀=-2param.L/param.N;
    (η2,v2,mesh2)=SolitaryWaveWhithamBoussinesq(param2;x₀=x₀)

    # Build in a second way the solitary wave
    init = SolitaryWB(param2)
    init2 = SolitaryWB(param2.c,ϵ=param2.ϵ,μ=param2.μ,N=2^10)

    # Check they correspond
    @test init.η(mesh.x)≈init2.η(mesh.x)
    @test init.v(mesh.x)≈init2.v(mesh.x)
    @test init.η(mesh.x)≈η
    @test init.v(mesh.x)≈v
    @test init.η(mesh2.x .-x₀)≈η2
    @test init.v(mesh2.x .-x₀)≈v2
    @test η2[end÷4:3*end÷4-1]≈η[1:end]
    @test v2[end÷4:3*end÷4-1]≈v[1:end]

    # Solve the corresponding initial-value problem
    param=merge(param2,(T=1e-2,dt=1e-3))
    model = WhithamBoussinesq(param)
    problem = Problem(model, init, param)
    solve!(problem;verbose=false)
    # Check the solution is translated
    (η,v,x,t)=solution(problem)
    @test η≈init.η(x.-param.c*param.T)
    @test v≈init.v(x.-param.c*param.T)

    # Do it again with other Whitham-Boussinesq model
    init = SolitaryWB(param;α=1/2)
    model = WhithamBoussinesq(param;α=1/2)
    problem = Problem(model, init, param)
    solve!(problem;verbose=false)
    (η,v,x,t)=solution(problem)
    @test η≈init.η(x.-param.c*param.T)
    @test v≈init.v(x.-param.c*param.T)

    # Do it again with regular Boussinesq
    init = SolitaryWB(param;Boussinesq=true,a=-1/4,b=1/2)
    model = WhithamBoussinesq(param;Boussinesq=true,a=-1/4,b=1/2)
    problem = Problem(model, init, param)
    solve!(problem;verbose=false)
    (η,v,x,t)=solution(problem)
    @test η≈init.η(x.-param.c*param.T)
    @test v≈init.v(x.-param.c*param.T)


end

#--- tests on Whitham-Green-Naghdi solitary wave
@testset "Whitham-Green-Naghdi solitary wave" begin
    # Build a solitary wave with different parameters
    param = (ϵ  = 0.1, μ = 0.1, c=1.1, L = 10, N=2^8)
    waves = []
    push!(waves,SolitaryWaveWhithamGreenNaghdi(param))
    push!(waves,SolitaryWaveWhithamGreenNaghdi(param;iterative=true))
    push!(waves,SolitaryWaveWhithamGreenNaghdi(param;
            iterative=true,max_iter = 15, tol = 1e-9, gtol = 1e-9, ktol = 1e-9))
    push!(waves,SolitaryWaveWhithamGreenNaghdi(param;dealias=1))
    push!(waves,SolitaryWaveWhithamGreenNaghdi(param;q=0.9,α=1/2,max_iter=30))
    push!(waves,SolitaryWaveWhithamGreenNaghdi(param;α=1/2))
    push!(waves,SolitaryWaveWhithamGreenNaghdi(param;method=1))
    push!(waves,SolitaryWaveWhithamGreenNaghdi(param;method=2))
    push!(waves,SolitaryWaveWhithamGreenNaghdi(param;method=3))
    push!(waves,SolitaryWaveWhithamGreenNaghdi(param;method=4))

    # Check they are consistent
    (η,u,v,mesh)=waves[1]
    for wave in waves
        (η2,u2,v2,mesh2)=wave
        @test η2≈η
        @test v2≈v
        @test u2≈u
    end

    # Build the translated solitary wave on a different mesh
    param2 = (ϵ  = 0.1, μ = 0.1, c=1.1, L = 20, N=2^9);x₀=-2param.L/param.N;
    (η2,u2,v2,mesh2)=SolitaryWaveWhithamGreenNaghdi(param2;x₀=x₀)

    # Build in a second way the solitary wave
    init = SolitaryWGN(param2)
    init2 = SolitaryWGN(param2.c,ϵ=param2.ϵ,μ=param2.μ,N=2^10)

    # Check they correspond
    @test init.η(mesh.x)≈init2.η(mesh.x)
    @test init.v(mesh.x)≈init2.v(mesh.x)
    @test init.η(mesh.x)≈η
    @test init.v(mesh.x)≈v
    @test init.η(mesh2.x .-x₀)≈η2
    @test init.v(mesh2.x .-x₀)≈v2
    @test η2[end÷4:3*end÷4-1]≈η[1:end]
    @test u2[end÷4:3*end÷4-1]≈u[1:end]
    @test v2[end÷4:3*end÷4-1]≈v[1:end]

    # Solve the corresponding initial-value problem
    param=merge(param2,(T=1e-2,dt=1e-3))
    model = WhithamGreenNaghdi(param)
    problem = Problem(model, init, param)
    solve!(problem;verbose=false)
    # Check the solution is translated
    (η,v,x,t)=solution(problem)
    @test η≈init.η(x.-param.c*param.T)
    @test v≈init.v(x.-param.c*param.T)

    # Do it again with regular Green-Naghdi
    init = SolitaryWGN(param;SGN=true)
    init2 = SolitarySGN(param)
    @test init.η(x)≈init2.η(x)
    @test init.v(x)≈init2.v(x)
    model = WhithamGreenNaghdi(param;SGN=true)
    problem = Problem(model, init, param)
    solve!(problem;verbose=false)
    (η,v,x,t)=solution(problem)
    @test η≈init.η(x.-param.c*param.T)
    @test v≈init.v(x.-param.c*param.T)


end

#--- tests on Whitham solitary wave
@testset "Whitham solitary wave" begin
    # Build a solitary wave with different parameters
    param = (ϵ  = 0.1, μ = 0.1, c=1.1, L = 10, N=2^9)
    waves = []
    push!(waves,SolitaryWaveWhitham(param;α=1))
    push!(waves,SolitaryWaveWhitham(param;iterative=true))
    push!(waves,SolitaryWaveWhitham(param;
            iterative=true,max_iter = 15, tol = 1e-9, gtol = 1e-9, ktol = 1e-9))
    push!(waves,SolitaryWaveWhitham(param;dealias=1))
    push!(waves,SolitaryWaveWhitham(param;q=0.9,α=1/2,max_iter=30))
    push!(waves,SolitaryWaveWhitham(param;α=1/2))

    # Check they are consistent
    (u,mesh)=waves[1]
    for wave in waves
        (u2,mesh2)=wave
        @test u2≈u
    end

    # Build the translated solitary wave on a different mesh
    param2 = (ϵ  = 0.1, μ = 0.1, c=1.1, L = 20, N=2^10);x₀=-2param.L/param.N;
    (u2,mesh2)=SolitaryWaveWhitham(param2;x₀=x₀,α=1)

    # Build in a second way the solitary wave
    init = SolitaryWhitham(param2)
    init2 = SolitaryWhitham(param2.c,ϵ=param2.ϵ,μ=param2.μ,N=2^10)

    # Check they correspond
    @test init.η(mesh.x)≈init2.η(mesh.x)
    @test init.v(mesh.x)≈init2.v(mesh.x)
    @test init.v(mesh.x)≈u
    @test init.v(mesh2.x .-x₀)≈u2
    @test u2[end÷4:3*end÷4-1]≈u[1:end]


end

#--- tests on randomly generated wave
@testset "Random initial data" begin
    param=(N=20,L=2);mesh=Mesh(param);x=mesh.x;
    # Build initial data
    inits=[]
    push!(inits,Random(param;L=1,s=Inf,λ=Inf,a=(1,1)))
    push!(inits,Random(param;L=1,s=2,λ=Inf,a=(1,1)))
    push!(inits,Random(param;L=1,s=Inf,λ=1,a=(1,1)))
    push!(inits,Random(param;L=1,s=Inf,λ=Inf,a=(-1,0)))


    # Check the data have correct size and none of them yields NaN
    for init in inits
        @test !any(isnan,init.η(x))
        @test !any(isnan,init.v(x))
        @test length(init.η(x))==param.N
        @test length(init.v(x))==param.N
    end
end
