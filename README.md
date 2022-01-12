# WaterWaves1D.jl

[![Build Status](https://github.com/WaterWavesModels/WaterWaves1D.jl/workflows/CI/badge.svg)](https://github.com/WaterWavesModels/WaterWaves1D.jl/actions)
<!-- [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://waterwavesmodels.github.io/WaterWaves1D.jl/stable) -->
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://waterwavesmodels.github.io/WaterWaves1D.jl/dev)

[![codecov](https://codecov.io/gh/WaterWavesModels/WaterWaves1D.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/WaterWavesModels/WaterWaves1D.jl)

## Installation

~~~
(v1.0) pkg> add https://github.com/WaterWavesModels/WaterWaves1D.jl.git
using WaterWaves1D
~~~

## Overview

`WaterWaves1D` provides a framework to compare several models for the propagation of unidimensional surface gravity waves.

Several models are already implemented, included (but not limited to) the water waves system, its truncated spectral expansion, the Green-Naghdi system, the Matsuno system...

An example of a possible usage of the code can be found below. More examples are available at the [examples](examples/) and [notebooks](notebooks/) repertories.


## Example

Define parameters of your problem
~~~
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
~~~

Define initial data
~~~
z(x) = exp.(-abs.(x).^4); # surface deformation
v(x) = 0*exp.(-x.^2);     # zero initial velocity
init = Init(z,v);         # generate the initial data with correct type
~~~

Set up initial-value problems for different models to compare
~~~
model1=WaterWaves(param,verbose=false) # The water waves system
model2=WWn(param;n=2,dealias=1,δ=1/10,verbose=false) # The quadratic model (WW2)
# type `?WaterWaves` or `?WWn` to see details and signification of arguments
problem1=Problem(model1, init, param, solver=RK4(model1)) ;
problem2=Problem(model2, init, param, solver=RK4(model2)) ;
~~~

Solve numerical time integration
~~~
solve!(problem1);
solve!(problem2);
~~~

Plot solutions at final time
~~~
plot_solution(problems;fourier=false)
~~~
![](./notebooks/Example.pdf)

Generate animation
~~~
anim = create_animation(problems;fourier=false,ylims=(-0.25,1))
import Plots.gif
gif(anim, "Example.gif", fps=15)
~~~
![](./notebooks/Example.gif)
