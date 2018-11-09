using Test
using FFTW
using LinearAlgebra
using ProgressMeter
using BenchmarkTools
using DeepWaterModels

epsilon = 1/2
N       = 2^10
L       = 10
T       = 5
dt      = 0.01

@show epsilon,N,L,T,dt

mesh   = Mesh(-L, L, N)
times  = Times(dt, T)
solver = RK4( N )

h  = zeros(Complex{Float64}, N)
u  = zeros(Complex{Float64}, N)

cheng = Cheng(mesh, epsilon)
matsuno = Matsuno(mesh, epsilon)

abstract type Model end

@testset "Parameters" begin

    param = Parameters( ϵ  = 1/2, 
                        N  = 2^12,
                        L  = 10,
                        T  = 5,
                        dt = 0.001)
    
    @test param.ϵ  == 0.5
    @test param.N  == 4096
    @test param.L  == 10
    @test param.T  == 5
    @test param.dt == 0.001

end

abstract type InitialData end

struct Bump <: InitialData  

    h :: Array{Complex{Float64},1}
    u :: Array{Complex{Float64},1}

    function Bump(p :: Parameters) 

    	mesh  = Mesh(-p.L, p.L, p.N)
    	h = exp.(-(mesh.x).^2)
    	u = zeros(Complex{Float64}, mesh.N)
    	new(h,u)

    end
end

struct Problem

    model   :: AbstractModel
    initial :: InitialData
    param   :: Parameters
    solver  :: TimeSolver
    data    :: Vector{Tuple{Vector{Complex{Float64}},
			    Vector{Complex{Float64}}}}

    function Problem(model   :: AbstractModel,
         	     initial :: InitialData,
         	     param   :: Parameters,
         	     solver  :: TimeSolver)

         data = [] 

         new(model,initial,param,solver,data)

    end
end

@testset "Initial data" begin

    param = Parameters( ϵ  = 1/2, 
                        N  = 2^12,
                        L  = 10,
                        T  = 5,
                        dt = 0.001)

    bump = Bump( param )
    @test true

end

@testset "Problem" begin

    param = Parameters( ϵ  = 1/2, 
                        N  = 2^12,
                        L  = 10,
                        T  = 5,
                        dt = 0.001)

    bump    = Bump( param )
    mesh    = Mesh(-param.L, param.L, param.N)
    cheng   = Cheng(mesh, param.ϵ)
    times   = Times(param.dt, param.T)
    solver  = RK4( param.N )
    problem = Problem( cheng, bump, param, solver )

    @test true

end

function solve!(model::AbstractModel, h, u, times::Times, solver::TimeSolver)
                
    prog = Progress(times.Nt,1) 
    
    model.data = []
  
    push!(model.data,(h,u))
    for l in range(1,times.Nt-1)
        
        dt = times.t[l+1]-times.t[l]
        
        step!( solver, model, h, u, dt)
    
        push!(model.data,(h,u))   
        next!(prog)
    end
            
end

h .= exp.(-mesh.x.^2)
u .= 0.0

h .= cheng.Pi .* fft(h)
u .= cheng.Pi .* fft(u)

solve!(cheng, h, u, times, solver )

@test !any(isnan,h)
@test !any(isnan,u)

h .= exp.(-mesh.x.^2)
u .= 0.0

h .= matsuno.Pi .* fft(h)
u .= matsuno.Pi .* fft(u)

solve!(matsuno, h, u, times, solver )

@test !any(isnan,h)
@test !any(isnan,u)
