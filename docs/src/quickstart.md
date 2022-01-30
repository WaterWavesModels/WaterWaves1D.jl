# Quickstart

# Notre premier exemple
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
