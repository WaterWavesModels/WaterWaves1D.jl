# Quickstart

`WaterWaves1D` provides a framework to study and compare several models for the propagation of unidimensional surface gravity waves (a.k.a. "water waves").

In the following examples, we
1. build a solitary wave profile for a fully dispersive "Whitham-Boussinesq" system;
2. check the numerical accuracy of the produced solitary wave and time integration solver;
3. integrate numerically the complete water waves system with the solitary wave profile as initial data.


## Build the solitary wave profile
```@example 1
using WaterWaves1D

param = ( ϵ  = 1/2,
          N  = 2^12,
          L  = 10.,
          T  = 5.,
          dt = 0.001,
          θ = 2.5 )

init    = BellCurve(param)
model   = Matsuno(param)
problem = Problem(model, init, param)

solve!( problem )

plot_solution(problem)

```
```@example 2
using WaterWaves1D
```
```@example 1

plot_solution(problem)
```

## Solve the Whitham-Boussinesq system

## Solve the water waves systen
