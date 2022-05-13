# WaterWaves1D.jl documentation



## Overview

[`WaterWaves1D.jl`](https://github.com/WaterWavesModels/WaterWaves1D.jl/) is a [Julia](https://julialang.org/) package providing a framework to study and compare several models for the propagation of unidimensional surface gravity waves (a.k.a. ["water waves"](background.md#Water-waves)).

Several models are already implemented, including ([but not limited to](background.md#Models)) the so-called water waves system, its truncated spectral expansion, the Green-Naghdi system, the Matsuno system, and so on. You may easily add your favorite one to the gang: see the [how-to guide](how-to.md#build-your-model).

## Installation

Enter the Pkg REPL by pressing `]` from the Julia REPL.
~~~
(v1.0) pkg> add https://github.com/WaterWavesModels/WaterWaves1D.jl.git
~~~

Once installed, load the package with
```@repl
using WaterWaves1D
```

## Quick start

A simple example is [documented here](example.md). More advanced examples can be found in the package's [examples](https://github.com/WaterWavesModels/WaterWaves1D.jl/tree/master/examples) and [notebooks](https://github.com/WaterWavesModels/WaterWaves1D.jl/tree/master/notebooks) folders.


## Developers

[`WaterWaves1D.jl`](https://github.com/WaterWavesModels/WaterWaves1D.jl/) is being developed by [Vincent Duchêne](https://perso.univ-rennes1.fr/vincent.duchene/) and [Pierre Navaro](https://github.com/pnavaro).

## Documentation contents

```@contents
Depth = 3
Pages = [
        "index.md",
        "background.md",
        "problems.md",
        "how-to.md",
        "example.md",
        "library.md"
        ]
```
