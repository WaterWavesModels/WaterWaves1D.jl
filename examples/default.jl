export plot_sol, plot_diff, anim
using WaterWaves1D, Plots

##########################################################################################
### Comparison of the water waves equations and quadratic WW2 model on a heap of water ###
##########################################################################################

#####################################################################################
# Gather parameters of the problem.
param = (
    # Physical parameters. Variables are non-dimensionalized as in Lannes, The water waves problem, isbn:978-0-8218-9470-5
    μ = 1,     # shallow-water dimensionless parameter
    ϵ = 1 / 4,   # nonlinearity dimensionless parameter
    # Numerical parameters
    N = 2^11,  # number of collocation points
    L = 10,    # half-length of the numerical tank (-L,L)
    T = 5,     # final time of computation
    dt = 0.01,  # timestep
);

#####################################################################################
# Define initial data (here, a "heap of water").

z(x) = exp.(-abs.(x) .^ 4); # surface deformation
v(x) = 0 * exp.(-x .^ 2);     # zero initial velocity
init = Init(z, v);         # generate the initial data with correct type

#####################################################################################
# Set up initial-value problems for different models to compare.

# Build models
model_WW = WaterWaves(param, verbose = false) # The water waves system
model_WW2 = WWn(param; n = 2, dealias = 1, δ = 1 / 10) # The quadratic model (WW2)
# Build problems
problem_WW = Problem(model_WW, init, param) ;
problem_WW2 = Problem(model_WW2, init, param) ;

#####################################################################################
# Integrate in time the initial-value problems.
solve!([problem_WW problem_WW2]);

# Conservation of mass, momentum, energy
Δ = [mass_diff(problem_WW2), momentum_diff(problem_WW2), energy_diff(problem_WW2)]

@info("Problems solved. 
Error in the conservation of mass : $(Δ[1]).
Error in the conservation of momentum : $(Δ[2]).
Error in the conservation of energy : $(Δ[3]).
Type `plot_sol()`, `plot_diff()` or `anim()` for illustrations.")

"""
    plot_sol()

Return the plot of the solutions of the two problems. 
See the documentation for advanced plotting recipes.

Requires the Plots package (`using Plots`).

# Example
```@example
using Plots
plt=plot_sol()
savefig(plt,"default_example.pdf")
```
"""
function plot_sol()
    return plot([problem_WW, problem_WW2])
end

"""
    plot_diff()

Return the plot of the difference between the solutions of the two problems. 
See the documentation for advanced plotting recipes.

Requires the Plots package (`using Plots`).

# Example
```@example
using Plots
plt=plot_sol()
savefig(plt,"default_example_difference.pdf")
```
"""
function plot_diff()
    return plot([(problem_WW, problem_WW2)], var = [:surface, :velocity])
end

"""
    anim(;name="default_example")

Generate an animation of the evolution of the solutions of the two problems. 

Requires the Plots package (`using Plots`).

# Example
```@example
using Plots
animation=anim()
gif(animation, "default_example.gif", fps=15)

```
"""
function anim()
    animation = @animate for t in LinRange(0, 5, 101)
        plot([problem_WW, problem_WW2]; T = t, ylims = (-0.5, 1))
    end
    return animation
end
nothing
