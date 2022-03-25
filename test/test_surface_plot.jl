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

model1 = WaterWaves(param) # The water waves system
model2 = WWn(param;n=2,dealias=1,δ=1/10) # The quadratic model (WW2)

problem1 = Problem(model1, init, param) ;
problem2 = Problem(model2, init, param) ;

solve!([problem1 problem2]; verbose=false);

plot(problem1)
plot!(problem2; legend = :bottomright)
