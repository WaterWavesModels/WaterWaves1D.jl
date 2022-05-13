using Test
import LinearAlgebra.norm
@testset "Deep water" begin
    η(x)=exp.(-x.^2); # Gaussian initial data for the surface deformation.
    v(x)=zero(x)      # we set the initial velocity as zero to avoid inconsistencies among different models.
    init=Init(η,v); 
    param_unstable = ( 
    # Physical parameters. Variables are non-dimensionalized as in Lannes, The water waves problem, isbn:978-0-8218-9470-5
    ϵ  = 1/4,    # nonlinearity dimensionless parameter
    μ  = Inf,    # inifinite-depth case
    # Numerical parameters
    N  = 2^11,   # number of collocation points
    L  = 2*π,    # half-length of the numerical tank (-L,L)
    T  = 1.5,      # final time of computation
    dt = 0.01,  # timestep
                );
    WW2_unstable = Problem( WWn(param_unstable,dealias = 1), init, param_unstable ) 
    solve!(WW2_unstable);
    WW2_rectified = Problem( WWn(param_unstable,dealias = 1,δ=0.01), init, param_unstable,label="rWW2" ) # δ is the strength of the regularization
    solve!(WW2_rectified);
    WW2η,=solution(WW2_unstable;T=1);rWW2η,=solution(WW2_rectified;T=1);
    errWW2=norm(WW2η-rWW2η)/sqrt(length(WW2η))
    @test errWW2 ≈ 4.21076939181685e-5
    WW2η,=solution(WW2_unstable)
    @test any(isnan,WW2η)

    param_heap = ( 
    # Physical parameters. Variables are non-dimensionalized as in Lannes, The water waves problem, isbn:978-0-8218-9470-5
    ϵ  = 1/4,    # nonlinearity dimensionless parameter
    μ  = Inf,    # infinite-depth case
    # Numerical parameters
    N  = 2^9,   # number of collocation points
    L  = 2*π,    # half-length of the numerical tank (-L,L)
    T  = 1,      # final time of computation
    dt = 0.01,  # timestep
    )
    init=Init(x->exp.(-x.^2),x->zero(x)); 
    WW_heap = Problem( WaterWaves(param_heap,dealias = 1, verbose=false), init, param_heap );
    rWW2_heap = Problem( WWn(param_heap,dealias = 1, δ=0.02), init, param_heap, label="rWW2" ) # δ is the strength of the regularization
    Matsuno_heap = Problem( Matsuno(param_heap,dealias = 1), init, param_heap )
    AkersNicholls_heap = Problem( AkersNicholls(param_heap,dealias = 1), init, param_heap )

    solve!(WW_heap);solve!(rWW2_heap);solve!(Matsuno_heap);solve!(AkersNicholls_heap);
    WWη,WWv,WWx,=solution(WW_heap)
    rWW2η,=solution(rWW2_heap,x=WWx);errWW2=norm(WWη-rWW2η)/sqrt(length(WWx))
    Matη,=solution(Matsuno_heap,x=WWx);errMat=norm(WWη-Matη)/sqrt(length(WWx))
    ANη,=solution(AkersNicholls_heap,x=WWx);errAN=norm(WWη-ANη)/sqrt(length(WWx))

    @test errWW2≈0.0008247817846683348
    @test errMat≈0.001381758315810658
    @test errAN ≈0.0007527210525751591

    param_rand = ( 
    # Physical parameters. Variables are non-dimensionalized as in Lannes, The water waves problem, isbn:978-0-8218-9470-5
    ϵ  = 1/4,    # nonlinearity dimensionless parameter
    μ  = Inf,    # infinite-depth case
    # Numerical parameters
    N  = 2^10,   # number of collocation points
    L  = 2*π,    # half-length of the numerical tank (-L,L)
    T  = 1,      # final time of computation
    dt = 0.01,  # timestep
    )
    mesh=Mesh(param_rand);x=mesh.x;
    init= Random(x;L=2,a=(1,0)); # we set the initial velocity as zero to avoid possible inconsistencies among different models.
    WW_rand = Problem( WaterWaves(param_rand,dealias = 1, verbose=false), init, param_rand );
    rWW2_rand = Problem( WWn(param_rand,dealias = 1, δ=0.05), init, param_rand, label="rWW2" );
    Matsuno_rand = Problem( Matsuno(param_rand,dealias = 1), init, param_rand );
    AkersNicholls_rand = Problem( AkersNicholls(param_rand,dealias = 1), init, param_rand );

    solve!(WW_rand);solve!(rWW2_rand);solve!(Matsuno_rand);solve!(AkersNicholls_rand);

    WWη,WWv,WWx,=solution(WW_rand)
    rWW2η,=solution(rWW2_rand,x=WWx);errWW2=norm(WWη-rWW2η)/sqrt(length(WWx))
    Matη,=solution(Matsuno_rand,x=WWx);errMat=norm(WWη-Matη)/sqrt(length(WWx))
    ANη,=solution(AkersNicholls_rand,x=WWx);errAN=norm(WWη-ANη)/sqrt(length(WWx))

    @test 1e-4< errWW2 <2*1e-3
    @test 1e-4< errMat <2*1e-3 
    @test 1e-4< errAN  <2*1e-3  
end