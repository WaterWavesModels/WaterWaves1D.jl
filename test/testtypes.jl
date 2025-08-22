using Test

@testset "Parameters" begin
    param = ( ϵ  = 1/2, μ = 1)
    paramX= ( N  = 2^8, L  = 10)
    paramT= ( T  = 0.1, dt = 0.01)

    @test param.ϵ  == 0.5
    @test param.μ  == 1
    @test paramX.N  == 256
    @test paramX.L  == 10
    @test paramT.T  == 0.1
    @test paramT.dt == 0.01
end
#--- test meshes


function the_mesh(m::Mesh)
    [m.N, m.xmin, m.xmax, m.dx, m.x, m.kmin, m.kmax, m.dk, m.k]
end

@testset "Tests on Mesh" begin

    m1a = Mesh( (xmin=-2*π,xmax=π,N=3) )
    m1b = Mesh( [-2*π -π 0 ] )

    m2a = Mesh( (xmin=-2*π,xmax=2π,N=2) )
    m2b = Mesh( [-2*π 0 ] )

    @test the_mesh(m1a) ≈ [3,-2π,π,π,[-2π; -π; 0],-2/3,2/3,2/3,[0; 2/3; -2/3]]
    @test the_mesh(m1a) == the_mesh(m1b)

    @test the_mesh(m2a) ≈ [2,-2π,2π,2π,[-2π; 0],-1//2,0,1//2,[0; -1/2]]
    @test the_mesh(m2a) == the_mesh(m2b)

end

#--- test times

function the_times(t::Times)
    [t.Nc, t.Ns, t.ns, t.tfin, t.dt, t.tc, t.ts]
end

@testset "Tests on times" begin
    t1a = Times( 0:0.1:1; ns=3 )
    t1b = Times( (dt=0.1 , T=1); ns=3 )

    t2a = Times( 0:1//10:1 ; Ns=3 )
    t2b = Times( (dt=0.1 , T=1); Ns=3 )

    @test the_times(t1a) == [11,4,[3,3,3],1.,0.1,0:0.1:1,[0.0, 0.3, 0.6, 0.9]]
    @test the_times(t1a) == the_times(t1b)

    @test the_times(t2a) == [11,4,[3,4,3],1,0.1,0:0.1:1,[0.0, 0.3, 0.7, 1.0]]
    @test the_times(t2a) == the_times(t2b)
end

#--- test data
@testset "Tests on data" begin
    v1 = [1; 3; 5];v2 = [2; 4; 6];
    M = [ v1 v2 ];
    data=Data(M)

    @test data.U==[[1 2 ; 3 4 ; 5 6 ]]
    @test data.datalength==length(v1)
    @test data.datasize==2
end

#--- test problems
@testset "Tests on problems" begin
    param = (L=10, N=20, ϵ = 0.1, μ = 0.1, T=1, dt=1/50)
    init = Random(param)
    solver=Euler_naive()
    # solver=RK4(model1)  # the tests will fail with this solver, since it uses pre-allocation
    # the tests also fail if we define  model=WWn(param;n=1), and use model to build the problems
    pb1  = Problem(WWn(param;n=1,mesh=Mesh(param)), init, param; solver=solver)
    pb1a = Problem(WWn(param;n=1), init, param; solver=solver)
    pb1b = Problem(WWn(param;n=1), init, param; solver=Euler_naive())
    pb2  = Problem(WWn(param;n=2), init, merge(param,(ns=2,)); solver=solver)
    pb2a = Problem(WWn(param;n=2), init, merge(param,(ns=10,)); solver=solver)
    pb2b = Problem(WWn(param;n=2), init, merge(param,(Ns=30,)); solver=solver)

    solve!([pb1a pb2a pb1b pb2b];verbose=false);
    solve!(pb1);solve!(pb2;verbose=true);

    for i in 1:2
        @test solution(pb1)[i]==solution(pb1a)[i]
        @test solution(pb1)[i]==solution(pb1b)[i]
        @test solution(pb2)[i]==solution(pb2a)[i]
        @test solution(pb2)[i]==solution(pb2b)[i]
    end
end

#--- test init
@testset "Tests on init" begin
    param = (L=10, N=2^8)
    mesh=Mesh(param);x=mesh.x;
    f(x) = exp.(-x.^2);g(x) = x.*exp.(-(x.-1).^2);
    inits=InitialData[]
    push!(inits, Init(f , g) )
    push!(inits, Init( mesh, f(x) , g(x); fast = false ) )
    push!(inits, Init( x, f(x) , g(x); fast = false ) )
    init = inits[1]
    for i in 2:3
        @test init.η([π])≈inits[i].η([π])
        @test init.v([π])≈inits[i].v([π])
    end
end
