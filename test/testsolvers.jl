using Test

#--- parameters
para = ( ϵ  = 0.1, μ = 0.1)  # physical parameters
paraX= ( N  = 4, L  = π)   # mesh with 6 collocation points: x=[-π, -π/2, 0.0, π/2]
paraT= ( T  = 1e-2, dt = 1e-3) # timegrid with ten instants: t=[0.0:1.0:10.0]/10^3
param = merge(para,paraX)  # used to construct models
parap = merge(paraX,paraT) # used to construct problems

#--- initial data
init     = Init(x->cos.(x),x-> sin.(x) )

#--- model
model=WWn(param)

#--- reference problem
pb0 = Problem( model, init, parap )
solve!(pb0,verbose=false)

#--- tests on Runge-Kutta 4 solvers
@testset "RK4 solvers" begin
    # different ways to build RK4 solver
    solvers = []
    push!(solvers, RK4(model) )
    push!(solvers, RK4(paraX) )
    push!(solvers, RK4(paraX.N) )
    push!(solvers, RK4(model.mapto(init)) )
    push!(solvers, RK4((paraX.N,2)) )
    push!(solvers, RK4_naive() )

    # check all RK4 solvers generate the same data
    for solver in solvers
        pb = Problem( model, init, parap;
                    solver = solver );
        solve!(pb,verbose=false)
        @test pb.data.U == pb0.data.U
    end
end

#--- tests on explicit Euler solvers
@testset "explicit Euler solver" begin
    # build explicit Euler solver
    pb1 = Problem( model, init, parap;
                solver = Euler(model) )
    solve!(pb1,verbose=false)
    # check explicit Euler is of order ≈ dt*T
    order = paraT.dt*paraT.T
    @test isapprox(pb0.data.U[end] , pb1.data.U[end] , rtol = order)
    @test !isapprox(pb0.data.U[end] , pb1.data.U[end] , rtol = order/4)


    # different ways to build explicit Euler solver
    solvers = []
    push!(solvers, Euler(paraX) )
    push!(solvers, Euler(paraX.N) )
    push!(solvers, Euler(model.mapto(init)) )
    push!(solvers, Euler((paraX.N,2)) )
    push!(solvers, Euler_naive() )

    # check all Euler solvers generate the same data
    for solver in solvers
        pb = Problem( model, init, parap;
                    solver = solver );
        solve!(pb,verbose=false)
        @test pb.data.U == pb1.data.U
    end
end

#--- tests on symplectic Euler solvers
@testset "symplectic Euler solver" begin
    # build symplectic Euler solver
    pb1 = Problem( model, init, parap;
                solver = EulerSymp(model) )
    solve!(pb1,verbose=false)
    # check symplectic Euler is of order ≈ dt*T
    order = paraT.dt*paraT.T
    @test isapprox(pb0.data.U[end] , pb1.data.U[end] , rtol = order)
    @test !isapprox(pb0.data.U[end] , pb1.data.U[end] , rtol = order/4)


    # different ways to build symplectic Euler solver
    solvers = []
    push!(solvers, EulerSymp(paraX) )
    push!(solvers, EulerSymp(paraX.N ) )
    push!(solvers, EulerSymp(model.mapto(init)) )

    # check all symplectic Euler solvers generate the same data
    for solver in solvers
        pb = Problem( model, init, parap;
                    solver = solver );
        solve!(pb,verbose=false)
        @test pb.data.U == pb1.data.U
    end

    # change number of iterations in the implicit step
    N=5
    pb2 = Problem( model, init, parap;
                solver = EulerSymp(model, Niter=N) )
    solve!(pb2,verbose=false)

    # check the difference is of order ≈ (dt)^(N+1)*T
    order = (paraT.dt)^(N+1)*paraT.T
    @test isapprox(pb2.data.U[end] , pb1.data.U[end], rtol = order )

    # change the equations solved by implicit method
    N=5
    pb3 = Problem( model, init, parap;
                solver = EulerSymp(model, implicit=2) )
    solve!(pb3,verbose=false)

    # check the difference is of order ≈ dt*T
    order = paraT.dt*paraT.T
    @test isapprox(pb3.data.U[end] , pb1.data.U[end], rtol = order )
    @test !isapprox(pb0.data.U[end] , pb1.data.U[end] , rtol = order/4)
end
