using Test
import LinearAlgebra.norm
@testset "Notebook: deep water" begin
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
    WW2η,=solution(WW2_unstable,T=1);rWW2η,=solution(WW2_rectified,T=1);
    errWW2=norm(WW2η-rWW2η)/sqrt(length(WW2η))
    @test errWW2 < 1e-4 && errWW2 > 1e-5
    WW2η,=solution(WW2_unstable)
    @test any(isnan,WW2η) && !any(isnan,rWW2η)

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

@testset "Notebook: Example" begin
    param = ( 
    # Physical parameters. Variables are non-dimensionalized as in Lannes, The water waves problem, isbn:978-0-8218-9470-5
    μ  = 1,     # shallow-water dimensionless parameter
    ϵ  = 1/4,   # nonlinearity dimensionless parameter
    # Numerical parameters
    N  = 2^9,  # number of collocation points
    L  = 10,    # half-length of the numerical tank (-L,L)
    T  = 5,     # final time of computation
    dt = 0.01,  # timestep
                );

    z(x) = exp.(-x.^4); # surface deformation
    v(x) = 0*z(x);      # zero initial velocity 
    init = Init(z,v);   # generate the initial data with correct type

    model1=WaterWaves(param,verbose=false) # The water waves system
    model2=WWn(param;n=2,dealias=1,δ=1/10) # The quadratic model (WW2)
    # type `?WaterWaves` or `?WWn` to see details and signification of arguments
    problem1=Problem(model1, init, param, solver=RK4(model1));
    problem2=Problem(model2, init, param, solver=RK4(model2));

    solve!(problem1);solve!(problem2);

    s1,=solution(problem1);
    s2,=solution(problem2);
    err=norm(s1-s2)/sqrt(length(s1))

    @test err≈0.026407484058654316

end

@testset "Notebook: Full dispersion" begin

    param = ( 
    # Physical parameters. Variables are non-dimensionalized as in Lannes, The water waves problem, isbn:978-0-8218-9470-5
    μ = 1,     # shallow-water dimensionless parameter
    ϵ = 1,     # nonlinearity dimensionless parameter
    c = 1.05,  # velocity of the solitary wave
    # Numerical parameters
    N = 2^7,  # number of collocation points
    L = 75    # half-length of the numerical tank (-L,L)
                    );
    # Compute the WGN solitary wave with velocity c (type `?SolitaryWaveWhithamGreenNaghdi` for more details)
    (ηWGN,uWGN,vWGN,mesh) = SolitaryWaveWhithamGreenNaghdi(param,verbose=false);
    # Compute the SGN solitary wave with velocity c (there is an exact formula)
    (ηSGN,uSGN,vSGN) = SolitaryWaveSerreGreenNaghdi(param);

    @test norm(ηWGN-ηSGN)/sqrt(param.N) ≈ 0.00018818525806149597
    @test norm(uWGN-uSGN)/sqrt(param.N) ≈ 0.00018061723014529588
    @test norm(vWGN-vSGN)/sqrt(param.N) ≈ 0.00018747696929027137



    T=2*param.L/param.c # final time of computation
    dt=T/10^3   # timestep
    init = Init(mesh,ηWGN,vWGN);           # use solitary wave as initial data
    model = WhithamGreenNaghdi(param)      # define the model 
    pb = Problem(                          # set up the initial-value problem to be solved
        model,
        init, 
        merge(param,(T=T,dt=dt)), 
        solver=RK4(model)   ) 
    solve!(pb)       

    ηfin,=solution(pb)
    err=norm(ηWGN-ηfin)/sqrt(param.N)

    @test err ≈ 2.5500799396961685e-6

    T=2*param.L/param.c # final time of computation
    dt=T/10^3   # timestep
    initWGN = Init(mesh,ηWGN,vWGN);       # WGN solitary wave as initial data
    initSGN = Init(mesh,ηSGN,vSGN);       # SGN solitary wave as initial data
    model = WaterWaves(param, dealias=0, verbose=false) # The water waves system
    pbSGN = Problem(                         # initial-value problem with SGN initial data
        model,
        initSGN, 
        merge(param,(T=T,dt=dt)), 
        solver=RK4(model)   ) 
    pbWGN = Problem(                         # initial-value problem with WGN initial data
        model,
        initWGN, 
        merge(param,(T=T,dt=dt)), 
        solver=RK4(model)   ) 
    solve!(pbSGN);solve!(pbWGN); 

    ηSGNfin,vSGNfin,xSGNfin=solution(pbSGN)
    ηWGNfin,vWGNfin,xWGNfin=solution(pbWGN)
    ηSGNinit,vSGNinit,xSGNinit=solution(pbSGN,T=0)
    ηWGNinit,vWGNinit,xWGNinit=solution(pbWGN,T=0)

    errSGN = norm(ηSGNinit-ηSGNfin)/sqrt(param.N)
    errWGN = norm(ηWGNinit-ηWGNfin)/sqrt(param.N)

    @test errSGN ≈ 0.00043699621029233914
    @test errWGN ≈ 7.999694318031115e-5

end

@testset "Notebook: Shallow water" begin

    param = ( 
        # Physical parameters. Variables are non-dimensionalized as in Lannes, The water waves problem, isbn:978-0-8218-9470-5
        μ  = 0.1,   # shallow-water dimensionless parameter
        ϵ  = 0.25,   # nonlinearity dimensionless parameter
        # Numerical parameters
        N  = 2^8,  # number of collocation points
        L  = π,    # half-length of the numerical tank (-L,L)
        T  = 1,    # final time of computation
        dt = 0.01, # timestep
                    );

    init= Init(x->exp.(-x.^4),zero); # we set the initial velocity as zero to avoid inconsistencies among different models.

    WW = Problem( WaterWaves(param,dealias = 1, verbose=false), init, param ) 
    GN = Problem( SerreGreenNaghdi(param; dealias = 0, iterate = true, precond = true), init, param )
    NH = Problem( NonHydrostatic(param; dealias = 0, iterate = true, precond = true), init, param )
    SRD= Problem( SquareRootDepth(param; dealias = 0, iterate = true, precond = true), init, param )
    WGN = Problem( WhithamGreenNaghdi(param; dealias = 0, iterate = true, precond = true), init, param )
    IK = Problem( IsobeKakinuma(param; dealias = 0, iterate = true, precond = true), init, param )


    solve!(WW);solve!(GN);solve!(NH);solve!(SRD);solve!(WGN);solve!(IK);

    WWη,WWv,WWx,=solution(WW)
    GNη,=solution(GN,x=WWx);errGN=norm(WWη-GNη)/norm(WWη)
    SRDη,=solution(SRD,x=WWx);errSRD=norm(WWη-SRDη)/norm(WWη)
    NHη,=solution(NH,x=WWx);errNH=norm(WWη-NHη)/norm(WWη)
    WGNη,=solution(WGN,x=WWx);errWGN=norm(WWη-WGNη)/norm(WWη)
    IKη,=solution(IK,x=WWx);errIK=norm(WWη-IKη)/norm(WWη)



    @test  errGN  ≈ 0.019492085082947582
    @test  errSRD ≈ 0.032762377673748344
    @test  errNH  ≈ 0.019380927144301712
    @test  errWGN ≈ 0.004015850164345053   
    @test  errIK  ≈ 0.002501413432973788

end

@testset "Notebook: Hammack-Segur" begin
    d  = 0.1;  # depth of the basin (in m)
    g  = 9.81; # gravitational acceleration
    λ  = 0.1;  # typical horizontal length (=d to simplify)
    T  = 40;   # final time (in s)
    L  = 156.16;# half-length of the numerical tank (in m)  
    param = ( 
        # Physical parameters. Variables are non-dimensionalized as in Lannes, The water waves problem, isbn:978-0-8218-9470-5
        μ  = 1,    # shallow-water dimensionless parameter
        ϵ  = 1,    # nonlinearity dimensionless parameter
        # Numerical parameters
        N  = 2^12, # number of collocation points
        L  = L/d,       # half-length of the numerical tank ( = 156.16 m)
        T  = T*sqrt(g*d)/λ, # final time of computation (= 50 s)
        dt = T*sqrt(g*d)/λ/10^3, # timestep
                    )
    using Elliptic
    a = 0.005 # half-amplitude of the initial data (in m)
    sn(x) = getindex.(ellipj.(9.25434*x*d,0.9999^2),1);  # regularized step function
    η0(x) = (-1/2*a.+1/2*a*sn(x)).*(x*d.<.61).*(x*d.>-1.83)/d;
    init = Init(x->2*η0(x),x->zero(x));         # generate the initial data with correct type

    model_WW2=WWn(param;n=2,dealias=1,δ=1/10) # The quadratic water waves model (WW2)
    #model_SGN=SerreGreenNaghdi(param) # The Serre-Green-Naghdi model (SGN)
    #model_WGN=WhithamGreenNaghdi(param) # The fully dispersive Whitham-Green-Naghdi model (WGN)
    # type `?WaterWaves` or `?WWn`, etc. to see details and signification of arguments
    WW2=Problem(model_WW2, init, param) ;
    #SGN=Problem(model_SGN, init, param) ;
    #WGN=Problem(model_WGN, init, param) ;

    solve!(WW2);#solve!(SGN);solve!(WGN);

    using HDF5
    data2a = h5read("./HammackSegur.h5","data2a")
    data2b = h5read("./HammackSegur.h5","data2b")
    data2c = h5read("./HammackSegur.h5","data2c")
    data2d = h5read("./HammackSegur.h5","data2d")
    data2e = h5read("./HammackSegur.h5","data2e")
    

    function gauge(p::Problem;x=0,T=nothing)
        if T==nothing
            times=p.times.ts
        elseif T[1]==T
            times=[T]
        else
            times=T
        end
        times/sqrt(g*d)*λ,[solution(p,T=ti,x=[x])[1][1] for ti in times]*d
    end
    nothing

    mesh_of_times = Mesh(WW2.times.ts /sqrt(g*d)*λ)

    ta=data2a[:,1]/sqrt(g*d)*λ
    ηa=data2a[:,2]*d*2/3
    ga=interpolate(mesh_of_times,gauge(WW2,x=1)[2],ta)

    @test norm(ηa-ga)/norm(ηa) ≈ 0.1672677795663852

    tb=(data2b[:,1].+50)/sqrt(g*d)*λ
    ηb=data2b[:,2]*d*2/3
    gb=interpolate(mesh_of_times,gauge(WW2,x=51)[2],tb)

    @test norm(ηb-gb)/norm(ηb) ≈ 0.13526329598410575


    tc=(data2c[:,1].+100)/sqrt(g*d)*λ
    ηc=data2c[:,2]*d*2/3
    gc=interpolate(mesh_of_times,gauge(WW2,x=101)[2],tc)

    @test norm(ηc-gc)/norm(ηc) ≈ 0.2599790311975911


    td=(data2d[:,1].+150)/sqrt(g*d)*λ
    ηd=data2d[:,2]*d*2/3
    gd=interpolate(mesh_of_times,gauge(WW2,x=151)[2],td)

    @test norm(ηd-gd)/norm(ηd) ≈ 0.4981462329721109


    te=(data2e[:,1].+200)/sqrt(g*d)*λ
    ηe=data2e[:,2]*d*2/3
    ge=interpolate(mesh_of_times,gauge(WW2,x=201)[2],te)

    @test norm(ηe-ge)/norm(ηe) ≈ 0.4921927105753341

end
