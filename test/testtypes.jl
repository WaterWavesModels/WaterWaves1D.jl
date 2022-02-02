using Test

#--- test meshes
m1a = Mesh( -2*π , π , 3 )
m1b = Mesh( (xmin=-2*π,xmax=π,N=3) )
m1c = Mesh( [-2*π -π 0 ] )

m2a = Mesh( -2*π , 2*π , 2 )
m2b = Mesh( (xmin=-2*π,xmax=2π,N=2) )
m2c = Mesh( [-2*π 0 ] )
m2d = Mesh( 2π , 2 )


function the_mesh(m::Mesh)
    [m.N, m.xmin, m.xmax, m.dx, m.x, m.kmin, m.kmax, m.dk, m.k]
end

@testset "Tests on Mesh" begin
    @test the_mesh(m1a) ≈ [3,-2π,π,π,[-2π; -π; 0],-2/3,2/3,2/3,[0; 2/3; -2/3]]
    @test the_mesh(m1a) == the_mesh(m1b)
    @test the_mesh(m1a) == the_mesh(m1c)

    @test the_mesh(m2a) ≈ [2,-2π,2π,2π,[-2π; 0],-1//2,0,1//2,[0; -1/2]]
    @test the_mesh(m2a) == the_mesh(m2b)
    @test the_mesh(m2a) == the_mesh(m2c)
    @test the_mesh(m2a) == the_mesh(m2d)

end

#--- test times
t1a = Times( 1//10 , 1; ns=3 )
t1b = Times( (dt=0.1 , T=1); ns=3 )

t2a = Times( 1//10 , 1; Ns=3 )
t2b = Times( (dt=0.1 , T=1); Ns=3 )

function the_times(t::Times)
    [t.Nc, t.Ns, t.ns, t.tfin, t.dt, t.tc, t.ts]
end

@testset "Tests on times" begin
    @test the_times(t1a) == [11,4,[3,3,3],1.,0.1,0:0.1:1,[0.0, 0.3, 0.6, 0.9]]
    @test the_times(t1a) == the_times(t1b)

    @test the_times(t2a) == [11,4,[3,4,3],1,0.1,0:0.1:1,[0.0, 0.3, 0.7, 1.0]]
    @test the_times(t2a) == the_times(t2b)
end

#--- test data
v1 = [1; 3; 5];v2 = [2; 4; 6];
M = [ v1 v2 ];
data=Data(M)
@testset "Tests on data" begin
    @test data.U==[[1 2 ; 3 4 ; 5 6 ]]
    @test data.datalength==length(v1)
    @test data.datasize==2
end
