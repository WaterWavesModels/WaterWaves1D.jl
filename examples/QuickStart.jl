# # Quick Start Example
#
#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/QuickStart.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/QuickStart    .ipynb)
#
# In this example we shall observe the disintegration of a heap of water using the water-waves system as well as a second-order small-steepness model.
#
# More advanced examples can be found in the package's [examples](https://github.com/WaterWavesModels/WaterWaves1D.jl/tree/master/examples) and [notebooks](https://github.com/WaterWavesModels/WaterWaves1D.jl/tree/master/examples/notebooks) folders.  Examples are also available through the functions [`examples_dir`](@ref WaterWaves1D.examples_dir), [`get_examples`](@ref WaterWaves1D.get_examples) and [`default_example`](@ref WaterWaves1D.default_example).
#
# ## Set up the initial-value problem
#
# First we define parameters of our problem.

using WaterWaves1D

param = (
     ## Physical parameters. Variables are non-dimensionalized as in Lannes, 
     ## The water waves problem, isbn:978-0-8218-9470-5
     μ  = 1,     # shallow-water dimensionless parameter
     ϵ  = 1/4,   # nonlinearity dimensionless parameter
     ## Numerical parameters
     N  = 2^10,  # number of collocation points
     L  = 10,    # half-length of the numerical tank (-L,L)
     T  = 5,     # final time of computation
     dt = 0.01,  # timestep
                 )

# Now we define [initial data solver](library.md#Initial-data) (the "heap of water"). The function [`Init`](@ref WaterWaves1D.Init) may take either functions, or vectors (values at collocation points) as arguments.

z(x) = exp.(-abs.(x).^4); # surface deformation
v(x) = zero(x);     # zero initial velocity
init = Init(z,v);         # generate the initial data with correct type

# Then we build the different [models](library.md#Models) to compare (see [`WaterWaves`](@ref WaterWaves1D.WaterWaves) and [`WWn`](@ref WaterWaves1D.WWn)).

model_WW = WaterWaves(param) # The water waves system
model_WW2 = WWn(param;n=2,dealias=1,δ=1/10) # The quadratic model (WW2)

# Finally we set up initial-value problems. Optionally, one may specify a [time solver](library.md#Solvers) to [`Problem`](@ref WaterWaves1D.Problem) (by default the standard explicit fourth order Runge Kutta method is used).

problem_WW = Problem(model_WW, init, param)
problem_WW2 = Problem(model_WW2, init, param)

# ## Solve the initial-value problem

# Solving the initial-value problems is as easy as [`solve!`](@ref WaterWaves1D.solve!).

solve!([problem_WW problem_WW2];verbose=false);

# ## Evaluate the precision of the numerical scheme

# The numerical discretization of preserved quantities, provided by the functions [`mass`](@ref WaterWaves1D.mass), [`momentum`](@ref WaterWaves1D.momentum), [`energy`](@ref WaterWaves1D.energy) (and [`mass_diff`](@ref WaterWaves1D.mass_diff), [`momentum_diff`](@ref WaterWaves1D.momentum_diff), [`energy_diff`](@ref WaterWaves1D.energy_diff) for the difference between initial and final time) provide valuable insights at the precision of a computed numerical solution.

Δ = [ mass_diff(problem_WW2), momentum_diff(problem_WW2), energy_diff(problem_WW2) ]

# ## Generate graphics

# Plot solutions at final time.

using Plots
plot([problem_WW, problem_WW2])

# Generate an animation.

@gif for t in LinRange(0,param.T,100)
    plot([problem_WW, problem_WW2], T = t)
    ylims!(-0.5, 1)
end

# See the [plot recipes](plot_recipes.md) for many more plotting possibilities.
