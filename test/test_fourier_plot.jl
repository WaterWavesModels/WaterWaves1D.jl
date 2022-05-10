using Plots
using WaterWaves1D

η(x) = exp.(-x.^2)
v(x) = zero(x)   
init = Init(η,v)

param1 = ( ϵ = 1/4, μ = Inf, N = 2^8, L = 2*π, T = 5, dt = 0.001 )

problem1 = Problem( WWn(param1, dealias = 1), init, param1 ) 

param2 = ( ϵ = 1/4, μ = Inf, N = 2^9, L = 2*π, T = 5, dt = 0.001 )

problem2 = Problem( WWn(param2, dealias = 1), init, param2 ) 

solve!([problem1, problem2])

plot(problem1; var = :surface)
plot!(problem2; var = :fourier)

#plot([problem1, problem2]; var = :fourier, label = :none)

