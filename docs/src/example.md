# Example

In this example we shall observe the disintegration of a heap of water using the water-waves system as well as a second-order small-steepness model.

More advanced examples can be found in the package's [examples](https://github.com/WaterWavesModels/WaterWaves1D.jl/tree/main/examples) and [notebooks](https://github.com/WaterWavesModels/WaterWaves1D.jl/tree/main/notebooks) section.

## Set up the initial-value problem

First we define parameters of our problem.

```@example 1
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
```

Now we define initial data (the "heap of water"). The function [`Init`](@ref WaterWaves1D.Init) may take either functions, or vectors (values at collocation points) as arguments.

```@example 1
z(x) = exp.(-abs.(x).^4); # surface deformation
v(x) = 0*exp.(-x.^2);     # zero initial velocity
init = Init(z,v);         # generate the initial data with correct type
```

Then we build the different models to compare (see [`WaterWaves`](@ref WaterWaves1D.WaterWaves) and [`WWn`](@ref WaterWaves1D.WWn)).

```@example 1
WW_model=WaterWaves(param) # The water waves system
WW2_model=WWn(param;n=2,dealias=1,δ=1/10) # The quadratic model (WW2)
```

Finally we set up initial-value problems. Optionally, one may specify a [time solver](index.md#Solvers) to [`Problem`](@ref WaterWaves1D.WaterWaves), by default the standard explicit fourth order Runge Kutta method is used.

```@example 1
WW_problem=Problem(WW_model, init, param) ;
WW2_problem=Problem(WW2_model, init, param) ;
```
## Solve the initial-value problem

Solving the initial-value problems is as easy as [`solve!`](@ref WaterWaves1D.solve!).

```@example 1
solve!([WW_problem WW2_problem];verbose=false);
```
## Generate graphics

Plot solutions at final time ([`plot_solution`](@ref WaterWaves1D.plot_solution) has many options).

```@example 1
plot_solution([WW_problem WW2_problem];fourier=false)
```

Generate an animation ([`plot_solution`](@ref WaterWaves1D.create_animation) has many options).


```@example 1
anim = create_animation([WW_problem WW2_problem];fourier=false,ylims=(-0.5,1))
gif(anim, "assets/example.gif", fps=15); nothing # hide
```

![](assets/example.gif)
