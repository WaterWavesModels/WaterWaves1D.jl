# Index

```@index
```


## Initial data

```@index
Pages   = [ "library/#initial-data.md" ]
```

```@autodocs
Modules = [WaterWaves1D]
#Filter = t -> typeof(t) === DataType && t <: InitialData # this filters out useful functions, so it is better to include specific pages
Pages = [
"CnoidalWaveSerreGreenNaghdi.jl",
"SolitaryWaveSerreGreenNaghdi.jl",
"SolitaryWaveWhithamGreenNaghdi.jl",
"SolitaryWaveWhithamBoussinesq.jl",
"SolitaryWaveWhitham.jl",
"Random.jl"]

```

## Models
```@index
Pages   = [ "/#models" ]
```


```@autodocs
Modules = [WaterWaves1D]
Filter = t -> typeof(t) === DataType && t <: AbstractModel
```

## Solvers
```@index
Pages   = [ "library#solvers.md" ]
```

```@autodocs
Modules = [WaterWaves1D]
Filter = t -> typeof(t) === DataType && t <: TimeSolver
```

## Structures

```@autodocs
Modules = [WaterWaves1D]
Pages   = [
"WaterWaves1D.jl",
"init.jl",
"mesh.jl",
"problem.jl",
"data.jl",
"times.jl"]
```


## Tools

```@autodocs
Modules = [WaterWaves1D]
Pages   = ["tools.jl"]
```

## Graphics

```@autodocs
Modules = [WaterWaves1D]
Pages   = ["figures.jl"]
```
