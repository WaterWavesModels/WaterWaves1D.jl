# Code basics

## Problems

A central object in [`WaterWaves1D.jl`](https://github.com/WaterWavesModels/WaterWaves1D.jl/) is the [`Problem`](@ref WaterWaves1D.Problem)) structure, which contains all information on a numerically discretized initial-value problem. In practice, a problem is generated as 
```
  problem = Problem( model, initial, times ; solver, label )
```
 where
 - `model` is related to the (spatially discretized) equation at stake. [Built-in models](index.md#Models) are typically generated as 
```
  model = TheModel( param ; kwargs )
```
where `param` is a `NamedTuple` containing relevant parameters of the model and of the spatial grid,  and `kwargs` are some optional arguments allowing some choices in the discretization (for instance [dealiasing](background.md#Pseudospectral-methods)), and a label for future references.
- `inital` is the couple of initial data. It can be generated for instance using the function [`Init`](@ref WaterWaves1D.Init) as  
```
  initial = Init( η, v )
```
where `η` and `v` are two functions (representing respectively the surface deformation and the derivative of the trave of the velocity potential at the surface). Alternatively, it can also be built from the values of these functions at equally-spaced collocation points.
- `times` contains relevant parameters of the time integration: in particular the final time  `T` and the time-step `dt`. It can be a `NamedTuple` with these informations or generated via the function [`Times`](@ref WaterWaves1D.Times).
- optionally, the time-solver `solver` can be provided (built-in solvers are the explicit Euler solver, [`Euler`](@ref WaterWaves1D.Euler) and [`Euler_naive`](@ref WaterWaves1D.Euler_naive), a symplectic Euler solver, [`EulerSymp`](@ref WaterWaves1D.EulerSymp), and the explicit Runge-Kutta 4 solver, [`RK4`](@ref WaterWaves1D.RK4) and [`RK4_naive`](@ref WaterWaves1D.RK4_naive)). By default the RK4 solver is used.
- optionally, a string `label` can be provided for future reference. It is inferred from the model if not provided.

The container `problem` then 
Once it has been built, `problem` contains, in addition to `model`, `initial`, `times`, `solver` and `label`, the raw data `data` (initially just the initial data). The initial-value problem is then numerically integrated (filling `data`) simply using
```
  solve!(problem)
```


## How to...

[`WaterWaves1D.jl`](https://github.com/WaterWavesModels/WaterWaves1D.jl/) is meant to be versatile,
and integrating new blocks to the package is easy. If you ever do so, please do not hesitate to contact the
[developers](home.md#Developers) to either get help, or report on your advances.

### build your model

Let us add the linear ([Airy](https://en.wikipedia.org/wiki/Airy_wave_theory)) water waves model whose equations are (using the notations introduced [here](background.md))
```math
  \left\{\begin{array}{l}
  ∂_tη+\frac{\tanh(\sqrt{μ} D)}{ν\sqrt{μ} D}∂_xv=0,\\[1ex]
  ∂_tv+∂_xη=0,
  \end{array}\right.
```
where we use the notation ``F(D)`` for the [action](https://en.wikipedia.org/wiki/Multiplier_(Fourier_analysis)) of pointwise multiplying by the function ``F`` in the Fourier space.

In a dedicated file we write
```julia
export Airy
mutable struct Airy <: AbstractModel
  label   :: String
  f!      :: Function
  mapto   :: Function
  mapfro  :: Function

  # We build our model here.

end
```

Our purpose is to provide
- `label`, a string used in subsequent informational messages, plots, etc.
- `f!`, a function to be called in explicit time-integration solvers such as [`Euler`](@ref WaterWaves1D.Euler) or [`RK4`](@ref WaterWaves1D.RK4) (one may provide other functions to be used with other solvers such as [`EulerSymp`](@ref WaterWaves1D.EulerSymp))
- `mapto`, a function which from  a couple `(η,v)` of type `InitialData` provides the raw data matrix on which computations are to be executed
- `mapfro`, a function which from raw data returns `(η,v,x)` where
    - `η` is the values of surface deformation at collocation points `x`;
    - `v` is the derivative of the trace of the velocity potential at `x`.


To this aim, in place of the commented line, we write
```julia
function Airy(param::NamedTuple; # param is a NamedTuple containing all necessary parameters
  label = "linear (Airy)"  # using a keyword argument allows the user to supersede the default label.
   )

  # Set up
  μ = param.μ
  if !in(:ν,keys(param)) # set default ν if it is not provided
    ν = min(1,1/√μ)
  else
    ν = param.ν
  end
  # Collocation points and Fourier modes
  m = Mesh(param)
  x, k = m.x, m.k
  # Fourier multipliers
  ∂ₓ	 = 1im * k            # Differentiation
  ∂ₓF₁ = 1im * tanh.(√μ*k)/(√μ*ν)
  # Pre-allocation
  fftη = zeros(Complex{Float64}, m.N)
  fftv = zeros(Complex{Float64}, m.N)

  # Evolution equations are ∂t U = f(U)
  function f!(U)
    fftη .= U[:,1]
    fftv .= U[:,2]

    U[:,1] .= -∂ₓF₁.*fftv
    U[:,2] .= -∂ₓ.*fftη
  end

  # Build raw data from physical data (discrete Fourier transform)
  function mapto(data::InitialData)
    U = [fft(data.η(x)) fft(data.v(x))]
  end

  # Return physical data `(η,v,x)` from raw data
  function mapfro(U)
    real(ifft(U[:,1])),real(ifft(U[:,2])),x
  end

  new( label, f!, mapto, mapfro )
end
```

A useful tool used in the above was the function [`Mesh`](@ref WaterWaves1D.Mesh).
It takes as argument a `NamedTuple` containing typically a number of points/modes
and the (half-)size of the domain and provides the vector of collocations points
and Fourier wavenumbers (see [this discussion](background.md#Pseudospectral-methods)).

The Airy model can now be built as follows
```julia
using WaterWaves1D
# include your file
model = Airy((μ=1,L=2π,N=2^8))
```

### build your initial data


### access to and manage your data

### plot your data

### save and load your data
