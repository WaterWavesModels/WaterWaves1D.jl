using Plots
using WaterWaves1D

param = (
    # Physical parameters. Variables are non-dimensionalized as in Lannes, The water waves problem, isbn:978-0-8218-9470-5
    μ  = 1,     # shallow-water dimensionless parameter
    ϵ  = 1/4,   # nonlinearity dimensionless parameter
    # Numerical parameters
    N  = 2^10,  # number of collocation points
    L  = 10,    # half-length of the numerical tank (-L,L)
    T  = 5,     # final time of computation
    dt = 0.01,  # timestep
                );
z(x) = exp.(-abs.(x).^4); # surface deformation
v(x) = zero(x);     # zero initial velocity
init = Init(z,v);         # generate the initial data with correct type

model = WaterWaves(param, tol = 1e-15) # The water waves system

problem = Problem(model, init, param)

solve!(problem; verbose=false)

plot(problem, var = [:surface, :velocity])
