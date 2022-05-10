using Plots
using WaterWaves1D

η(x) = exp.(-x.^2)
v(x) = zero(x)   
init = Init(η,v)

param = ( ϵ = 1/4, μ = Inf, N = 2^8, L = 2*π, T = 5, dt = 0.001 )

problem = Problem( WWn(param1, dealias = 1), init, param1 ) 

solve!(problem)

l = @layout [a ; b ; c]
p1 = plot(problem, var = :surface)
p2 = plot(problem, var = :velocity)
p3 = plot(problem, var = :fourier)
plot(p1, p2, p3, layout = l)
