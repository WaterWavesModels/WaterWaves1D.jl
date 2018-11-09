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

    model   :: Model
    initial :: InitialData
    param   :: Parameters
    solver  :: TimeSolver
    data    :: Vector{Tuple{Vector{Complex{Float64}},Vector{Complex{Float64}}}}

    function Problem(model   :: Model,
         	    initial :: InitialData,
         	    param   :: Parameters,
         	    solver  :: TimeSolver)

         data = [] 

         new(model,initial,param,solver,data)

    end
end

@testset "Parameters" begin

    params = Parameters( ϵ  = 1/2, 
                         N  = 2^12,
                         L  = 10,
                         T  = 5,
                         dt = 0.001)
    
    @test params.ϵ  == 0.5
    @test params.N  == 4096
    @test params.L  == 10
    @test params.T  == 5
    @test params.dt == 0.001

end

@testset "Initial data" begin

end

@testset "Problem" begin

    problem = Problem( cheng, initial, params, solver )

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


for model in [cheng, matsuno]
        
    h .= exp.(-mesh.x.^2)
    u .= 0.0

    h .= model.Pi .* fft(h)
    u .= model.Pi .* fft(u)
    
    solve!(model, h, u, times, solver )

end

@test true
