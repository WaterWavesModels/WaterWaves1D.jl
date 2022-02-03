# Quickstart

Define parameters of your problem

```@example 1
using Plots
using WaterWaves1D

param = (
    # Physical parameters. Variables are non-dimensionalized as in Lannes, The water waves problem, isbn:978-0-8218-9470-5
    μ  = 1,     # shallow-water dimensionless parameter
    ϵ  = 1/4,   # nonlinearity dimensionless parameter
    # Numerical parameters
    N  = 2^11,  # number of collocation points
    L  = 10,    # half-length of the numerical tank (-L,L)
    T  = 5,     # final time of computation
    dt = 0.01, # timestep
                );
```

Define initial data

```@example 1
z(x) = exp.(-abs.(x).^4); # surface deformation
v(x) = 0*exp.(-x.^2);     # zero initial velocity
init = Init(z,v);         # generate the initial data with correct type
```

Set up initial-value problems for different models to compare

```@example 1
model1=WaterWaves(param) # The water waves system
model2=WWn(param;n=2,dealias=1,δ=1/10) # The quadratic model (WW2)
# type `?WaterWaves` or `?WWn` to see details and signification of arguments
problem1=Problem(model1, init, param, solver=RK4(model1)) ;
problem2=Problem(model2, init, param, solver=RK4(model2)) ;
```

Solve numerical time integration

```@docs
solve!
```

```@example 1
solve!(problem1);
solve!(problem2);
```

Plot solutions at final time

```@example 1
plot_solution([problem1 problem2];fourier=false)
```

Generate an animation

```@example 1
anim = create_animation([problem1 problem2];fourier=false,ylims=(-0.5,1))
gif(anim, "assets/example.gif", fps=15) # hide
```

![](assets/example.gif)
