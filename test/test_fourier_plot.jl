using Plots
using WaterWaves1D

η(x) = exp.(-x.^2)
v(x) = zero(x)   
init = Init(η,v)

param = ( ϵ = 1/4, μ = Inf, N = 2^9, L = 2*π, T = 5, dt = 0.001 )
problem = Problem( WWn(param,dealias = 1), init, param ) 
solve!(problem)

plot(problem; fourier = true)

