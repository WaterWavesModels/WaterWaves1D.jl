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
    @test 4e-11 < a < 5e-11

    b=PlotSolitaryWaveWhitham(1.17)
    @test 0.0674 < b <  0.0675

end

@testset "Example: Solitary Wave Whitham-Boussinesq" begin
    include("../examples/SolitaryWaveWhithamBoussinesq.jl")

    b=IntegrateSolitaryWaveWhithamBoussinesq(N=2^6,dt=0.01/1.05)
    @test 9.9e-11 < b < 10e-11

    a=PlotSolitaryWaveWhithamBoussinesq(;c=1.1,α=1,L=20,N=2^9,μ=0.1,ϵ=0.1)
    @test 0.0316 < a < 0.0317

end

@testset "Example: Spectral Stability Green Naghdi" begin
    include("../examples/SpectralStabilityGreenNaghdi.jl")
    param = ell(0.3,0.005,0.75)
    σ,param,η,u,v,mesh=figspecCW(param;P=1,N=2^4)


    @test param == (h₀ = 0.3, h₁ = 0.3021875, h₂ = 0.305, a₀ = 0.3, a₁ = 0.0050000000000000044, H₀ = 0.3034497099913216, c = 0.166283361314354, λ = 5.189093852987355, m = 0.7499999999999991, κ = 0.3682704215595702)
    @test mesh == Mesh((N=2^4,L=param.λ))
    @test all((norm(σ),norm(η),norm(u),norm(v)) .≈ (92.75180153770384, 2.7862039840304336, 1.5268134085348848, 1.526814071924427))

    σ,η,u,v,mesh=figspecSW(1.2;L=10*π,N=2^6,M=10)
    @test mesh == Mesh((N=2^6,L=10*π))
    @test all((norm(σ),norm(η),norm(u),norm(v)) .≈ (95.64108213520367, 0.7411137440019041, 0.6632854772327135, 0.7064480663393101))

end

@testset "Study: Rectified WW2" begin
    include("../examples/StudyRectifiedWW2.jl")
    problem,blowup_time,blowup,error_energy=IntegrateWW2(init=1,μ=1,ϵ=0.1,L=20,N=2^6,T=1,dt = 0.01,dealias=0,δ=0,m=-1)
    @test (blowup_time,blowup) == (1.,false)
    @test 9e-7 < error_energy < 10e-7

    problem,blowup_time,blowup,error_energy=IntegrateWW2(init=1,μ=1,ϵ=10,L=20,N=2^14,T=0.1,dt = 0.01,dealias=0,δ=0,m=-1)
    @test blowup == true || blowup_time in (0.06,0.07)

end

@testset "Study: Whitham-Green-Naghdi" begin
    include("../examples/StudyWhithamGreenNaghdi.jl")

    norms = (13.564148081919193, 12.01160049359638)

    η,u,v,mesh = PlotSolitaryWaveWGN1(c=2,N=2^9,L=7*π,verbose=false,name=nothing) 
    @test all((norm(η),norm(v)) .≈ norms )
    @test mesh == Mesh((N=2^9,L=7*π))

    η,u,v,mesh = PlotSolitaryWaveWGN2(c=2,N=2^9,L=7*π,verbose=false,name=nothing) 
    @test all((norm(η),norm(v)) .≈ norms )
    @test mesh == Mesh((N=2^9,L=7*π))

    η,u,v,mesh = PlotSolitaryWaveWGN3(c=2,N=2^9,L=7*π,verbose=false,name=nothing) 
    @test all((norm(η),norm(v)) .≈ norms )
    @test mesh == Mesh((N=2^9,L=7*π))

    Jac,Jacstar,FFT,IFFT = PlotJacobianWGN(;c=2,L=10*π,N=2^6,SGN=false,verbose=false,name=nothing)
    @test norm(Jac) .≈ 98.38287473770893 || norm(Jacstar) .≈ 73.20464101226274
    @test norm(FFT) == 2^6 || norm(IFFT) == 1

    pb = IntegrateSolitaryWaveWGN(;SGN=true,c=2,N=2^8,L=8*π,T=1,dt=1/100,name=nothing)
    η0,u0,v0=SolitaryWaveSerreGreenNaghdi((c=2,N=2^8,L=8*π,ϵ=1,μ=1); x₀ = 2)
    η,v=solution(pb)
    @test norm(η-η0)/norm(η) < 2e-7 || norm(v-v0)/norm(v) < 4e-8

    pb = StabilitySolitaryWaveWGN(;p=2,c=2,N=2^6,L=5*π,T=1,dt=10/10^3,SGN=true,precond=true,iterate=true,dealias=0,name=nothing)
    η,v=solution(pb)
    @test norm(η) ≈ 5.720544559359206 || norm(v) ≈ 5.080826080626654

    pb = IntegrateWGN(2;δ=0.1,N=2^6,L=3*π,x₀=-3,T= 1,dt = 10/10^3,SGN=false,dealias=0,iterate=true,precond=true,name=nothing)
    η,v=solution(pb)
    @test norm(η) ≈ 1.0808748139030515 || norm(v) ≈ 0.9929455087249658
    

end

@testset "Study: SaintVenant" begin
    include("../examples/StudySaintVenant.jl")

    problem1 = IntegrateSV(;init=1,α=1.5,N=nothing,h₀=nothing,v₀=nothing,ϵ=1,L=π,M=2^6,T=0.05,dt =1e-5,dealias=true,smooth=false,Ns=1)
    η,v=solution(problem1)
    @test all((norm(η),norm(v)) .≈ (1.1384447011432766,0.13223692570972134))


    problem2 = IntegrateSV(;init=2,α=nothing,N=nothing,h₀=1/2,v₀=2,ϵ=1,L=π,M=2^6,T=0.05,dt =1e-5,dealias=true,smooth=true,Ns=1)
    η,v=solution(problem2)
    @test all((norm(η),norm(v)) .≈ (3.402895534611065, 11.157811912130184))

    problem3 = IntegrateSV(;init=3,α=nothing,N=nothing,h₀=nothing,v₀=nothing,ϵ=1,L=π,M=2^6,T=0.05,dt =1e-5,dealias=false,smooth=nothing,Ns=1)
    η,v=solution(problem3)
    @test all((norm(η),norm(v)) .≈ (1.26323755549213, 1.1565217449137426))

    problem2Da = IntegrateSV2D(; N=nothing, s=2, h₀=1/2,
                        u₀=1/2, v₀=-1/2, u₁=1, v₁=-1,
                        ϵ=1, L=π, M=2^5, T=0.01, dt=5e-5,
                        dealias=true, smooth=false,
                        hamiltonian=false, Ns=1,
                        label="non-hamiltonian")

    η,vx,vy=solution(problem2Da)
    @test all( (norm(η),norm(vx),norm(vy)) .≈ (7.999250008281116, 7.921618162378015, 8.081581836264434) )

    problem2Db = IntegrateSV2D(; N=nothing, s=2, h₀=1/2,
                        u₀=1/2, v₀=-1/2, u₁=1, v₁=-1,
                        ϵ=1, L=π, M=2^5, T=0.01, dt=5e-5,
                        dealias=true, smooth=false,
                        hamiltonian=true, Ns=1,
                        label="non-hamiltonian")

    η,vx,vy=solution(problem2Db)
    @test all( (norm(η),norm(vx),norm(vy)) .≈ (7.999225016708526, 7.92164586145311, 8.081605781638073) )


end