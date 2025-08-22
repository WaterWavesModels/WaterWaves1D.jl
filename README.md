# WaterWaves1D.jl

[![Build Status](https://github.com/WaterWavesModels/WaterWaves1D.jl/workflows/CI/badge.svg)](https://github.com/WaterWavesModels/WaterWaves1D.jl/actions)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://waterwavesmodels.github.io/WaterWaves1D.jl/stable/)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://waterwavesmodels.github.io/WaterWaves1D.jl/dev/)
[![codecov](https://codecov.io/gh/WaterWavesModels/WaterWaves1D.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/WaterWavesModels/WaterWaves1D.jl)
[![DOI](https://zenodo.org/badge/154723425.svg)](https://zenodo.org/badge/latestdoi/154723425)


## Installation

```julia
(v1.11) pkg> add WaterWaves1D
using WaterWaves1D
```

## Overview

`WaterWaves1D.jl` is a [Julia](https://julialang.org/) package providing a framework to study and compare several models for the propagation of unidimensional surface gravity waves (a.k.a. ["water waves"](https://waterwavesmodels.github.io/WaterWaves1D.jl/dev/background/#Water-waves)).

Several models are already implemented, including ([but not limited to](https://waterwavesmodels.github.io/WaterWaves1D.jl/dev/background/#Models)) the so-called water waves system, its truncated spectral expansion, the Green-Naghdi system, the Matsuno system, and so on. You may easily add your favorite one to the gang: see the [how-to guide](https://waterwavesmodels.github.io/WaterWaves1D.jl/dev/how-to/#build-your-model).

## Documentation

See [here](https://waterwavesmodels.github.io/WaterWaves1D.jl/dev/).


## Example

A simple example of a typical use of the package can be found below. More advanced examples are available at the [examples](examples/) and [notebooks](notebooks/) folders.



Gather parameters of the problem.
```julia
param = (
    # Physical parameters. Variables are non-dimensionalized as in Lannes, The water waves problem, isbn:978-0-8218-9470-5
    μ  = 1,     # shallow-water dimensionless parameter
    ϵ  = 1/4,   # nonlinearity dimensionless parameter
    # Numerical parameters
    N  = 2^11,  # number of collocation points
    L  = 10,    # half-length of the numerical tank (-L,L)
    T  = 5,     # final time of computation
    dt = 0.01,  # timestep
                );
```

Define initial data (here, a "heap of water").
```julia
z(x) = exp.(-abs.(x).^4); # surface deformation
v(x) = 0*exp.(-x.^2);     # zero initial velocity
init = Init(z,v);         # generate the initial data with correct type
```

Set up initial-value problems for different models to compare.
```julia
# Build models
model_WW=WaterWaves(param,verbose=false) # The water waves system
model_WW2=WWn(param;n=2,dealias=1,δ=1/10) # The quadratic model (WW2)
# Build problems
problem_WW=Problem(model_WW, init, param) ;
problem_WW2=Problem(model_WW2, init, param) ;
```

Integrate in time the initial-value problems.
```julia
solve!([problem_WW problem_WW2]);
```

Plot solutions at final time.
```julia
using Plots
plot([problem_WW, problem_WW2])
```
![](./notebooks/Example.png)

Generate an animation.
```julia
anim = @animate for t = LinRange(0,5,101)
    plot([problem_WW, problem_WW2];T=t,ylims=(-0.5,1))
end
gif(anim, "Example.gif", fps=15)
```
![](./notebooks/Example.gif)


## Developers

`WaterWaves1D.jl` is being developed by [Vincent Duchêne](https://perso.univ-rennes1.fr/vincent.duchene/) and [Pierre Navaro](https://github.com/pnavaro).

## Citing

The code is citable via [zenodo](https://zenodo.org). Please cite as:

> V. Duchêne, P. Navaro. WaterWaves1D.jl (Version v0.1.0). Zenodo.  [https://doi.org/10.5281/zenodo.7142921](https://doi.org/10.5281/zenodo.7142921)

## Publications

These preprints and publications use `WaterWaves1D.jl`.

1. Vincent Duchêne and Christian Klein, [*Numerical study of the Serre-Green-Naghdi equations and a fully dispersive counterpart*](https://doi.org/10.3934/dcdsb.2021300), Discrete Contin. Dyn. Syst. Ser. B, 27(10), pp. 5905–5933 (2022). DOI: [10.3934/dcdsb.2021300](https://doi.org/10.3934/dcdsb.2021300)
2. Vincent Duchêne and Benjamin Melinand, [*Rectification of a deep water model for surface gravity waves*](https://doi.org/10.2140/paa.2024.6.73), Pure Appl. Anal. 6, No. 1, pp. 73–128 (2024). DOI: [10.2140/paa.2024.6.73](https://doi.org/10.2140/paa.2024.6.73)
3. Arnaud Debussche, Etienne Mémin, and Antoine Moneyron, [*Derivation of Stochastic Models for
Coastal Waves*](https://doi.org/10.1007/978-3-031-70660-8_9), in B. Chapron et al. (eds.), Stochastic Transport in Upper Ocean Dynamics III,
Mathematics of Planet Earth 13 (2025). DOI:[10.1007/978-3-031-70660-8_9](https://doi.org/10.1007/978-3-031-70660-8_9)
4. Duchêne, Vincent; Marstrander, Johanna Ulvedal, [*The Fourier spectral approach to the spatial discretization of quasilinear hyperbolic systems*](https://doi.org/10.48550/arXiv.2507.00516), 
Preprint, arXiv:2507.00516 (2025). DOI: [10.48550/arXiv.2507.00516](https://doi.org/10.48550/arXiv.2507.00516) 

If you have a work that belongs to tis list, please [open a pull request](https://github.com/WaterWavesModels/WaterWaves1D.jl/pulls) to add it or let us know!

## Contributing

If you would like to contribute to the development of the package, please do not refrain from contacting us, for instance by [opening an issue](https://github.com/WaterWavesModels/WaterWaves1D.jl//issues/new).


## Related Julia packages

If this package does not fit your needs, you may consider

- [`SpectralWaves.jl`](https://github.com/mcpaprota/SpectralWaves.jl) Implementation of a spectral method for water waves *over topography*.
- [`GeophysicalFlows.jl`](https://github.com/FourierFlows/GeophysicalFlows.jl) Geophysical fluid dynamics leveraging the pseudo-spectral solver modules of the package [`FourierFlows.jl`](https://github.com/FourierFlows/FourierFlows.jl).
- [`DispersiveShallowWater.jl`](https://juliapackages.com/p/dispersiveshallowwater) Structure-preserving numerical methods for dispersive shallow water models.
- [`TrixiShallowWater.jl`](https://github.com/trixi-framework/TrixiShallowWater.jl) simulations of the shallow water equations with the discontinuous Galerkin method.
- [`Oceananigans.jl`](https://github.com/CliMA/Oceananigans.jl) Finite volume simulations of the nonhydrostatic and hydrostatic Boussinesq equations on CPUs and GPUs.