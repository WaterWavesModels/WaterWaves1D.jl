using Test
using FFTW
using LinearAlgebra
using ProgressMeter
using BenchmarkTools
using Plots
using DeepWaterModels

pyplot()

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

epsilon = 1/2
N       = 2^12
L       = 10
T       = 5
dt      = 0.001

@show epsilon,N,L,T,dt

mesh   = Mesh(-L, L, N)
times  = Times(dt, T)
solver = RK4( N )

h  = zeros(Complex{Float64}, N)
u  = zeros(Complex{Float64}, N)

models = [Cheng(mesh, epsilon), Matsuno(mesh, epsilon)]

for model in models
        
    h .= exp.(-mesh.x.^2)
    u .= 0.0

    h .= model.Pi .* fft(h)
    u .= model.Pi .* fft(u)
    
    solve!(model, h, u, times, solver )

end

fig(5, times, models, mesh)

savefig("test.png")

@test true
