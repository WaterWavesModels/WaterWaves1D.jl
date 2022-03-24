using Plots
using Printf
using WaterWaves1D
using FFTW

η(x)=exp.(-x.^2); # Gaussian initial data for the surface deformation.
v(x)=zero(x)      # we set the initial velocity as zero to avoid inconsistencies among different models.
init=Init(η,v); 

param = ( 
    # Physical parameters. Variables are non-dimensionalized as in Lannes, The water waves problem, isbn:978-0-8218-9470-5
    ϵ  = 1/4,    # nonlinearity dimensionless parameter
    μ  = Inf,    # inifinite-depth case
    # Numerical parameters
    N  = 2^9,    # number of collocation points
    L  = 2*π,    # half-length of the numerical tank (-L,L)
    T  = 5,      # final time of computation
    dt = 0.001,  # timestep
                );
problem = Problem( WWn(param,dealias = 1), init, param ) 
solve!(problem);

plot(problem, fourier = true)

