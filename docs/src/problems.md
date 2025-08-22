# Main architecture

A central object in [`WaterWaves1D.jl`](https://github.com/WaterWavesModels/WaterWaves1D.jl/) is the [`Problem`](@ref WaterWaves1D.Problem) structure, which contains all information on a numerically discretized initial-value problem. In practice, a problem is generated as 
```julia
problem = Problem( model, initial, times ; solver, label )
```
 where
 - `model` is related to the (spatially discretized) equation at stake. [Built-in models](library.md#Models) are typically generated as 
```julia
model = MyModel( param ; kwargs )
```
where `param` is a `NamedTuple` containing relevant parameters of the model and of the spatial grid,  and `kwargs` are some optional arguments allowing some choices in the discretization (for instance [dealiasing](background.md#Pseudospectral-methods)), and a label for future references.
- `inital` is the couple of initial data. It can be generated for instance using the function [`Init`](@ref WaterWaves1D.Init) as  
```julia
initial = Init( η, v )
```
where `η` and `v` are two functions (representing respectively the surface deformation and the derivative of the trace of the velocity potential at the surface). Alternatively, it can also be built from the values of these functions at equally-spaced collocation points.
- `times` contains relevant parameters of the time integration: in particular the final time  `T` and the time-step `dt`. It can be a `NamedTuple` with these informations or generated via the function [`Times`](@ref WaterWaves1D.Times).
- optionally, the time-solver `solver` can be provided (built-in solvers are the explicit Euler solver, [`Euler`](@ref WaterWaves1D.Euler) and [`Euler_naive`](@ref WaterWaves1D.Euler_naive), a symplectic Euler solver, [`EulerSymp`](@ref WaterWaves1D.EulerSymp), and the explicit Runge-Kutta 4 solver, [`RK4`](@ref WaterWaves1D.RK4) and [`RK4_naive`](@ref WaterWaves1D.RK4_naive)). By default the RK4 solver is used.
- optionally, a string `label` can be provided for future reference. It is inferred from the model if not provided.

Once it has been built, `problem` contains the raw data, `data` (initially just the initial data), in addition to `model`, `initial`, `times`, `solver` and `label`. The initial-value problem is then numerically integrated (filling `data`) with
```julia
solve!(problem)
```
