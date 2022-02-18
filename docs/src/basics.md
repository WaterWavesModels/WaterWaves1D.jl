# Code basics

Under construction

## How to...

[`WaterWaves1D.jl`](https://github.com/WaterWavesModels/WaterWaves1D.jl/) is meant to be versatile,
and integrating new blocks to the package is easy. If you ever do so, please do not hesitate to contact the
[developers](home.md#Developers) to either get help, or report on your advances.

### add your model

Let us add the linear ([Airy](https://en.wikipedia.org/wiki/Airy_wave_theory)) water waves model whose equations are (using the same notations as [here](background.md))
```math
  \left\{\begin{array}{l}
  ∂_tη+\frac{\tanh(\sqrt{μ} D)}{ν\sqrt{μ} D}∂_xv=0,\\[1ex]
  ∂_tv+∂_xη=0,
  \end{array}\right.
```
where ``F_1^μ=`` (here we use the notation ``F(D)`` for the [action](https://en.wikipedia.org/wiki/Multiplier_(Fourier_analysis)) of pointwise multiplying by the function ``F`` in the Fourier space).

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
- `f!` a function to be called in explicit time-integration solvers such as [`Euler`](@ref WaterWaves1D.Euler) or [`RK4`](@ref WaterWaves1D.RK4) (one may provide other functions to be used with other solvers such as [`EulerSymp`](@ref WaterWaves1D.EulerSymp))
- `WhithamBoussinesq.mapto` a function which from  a couple `(η,v)` of type `InitialData` provides the raw data matrix on which computations are to be executed
- `WhithamBoussinesq.mapfro` the inverse function which from raw data returns `(η,v)`.

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

  # Return physical data `(η,v)` from raw data
  function mapfro(U)
    real(ifft(U[:,1])),real(ifft(U[:,2]))
  end

  new( label, f!, mapto, mapfro )
end
```

A useful tool used in the above was the function [`Mesh`](@ref WaterWaves1D.Mesh).
It takes as argument a `NamedTuple` containing typically a number of points/modes
and the (half-)size of the domain and provides the vector of collocations points
and Fourier wavenumbers (see [this discussion](background.md#Pseudospectral-methods)).

The Airy model can now be built as follows
```@example 2
using WaterWaves1D
# include your file
model = Airy((μ=1,L=2π,N=2^8))
```

### add your initial data

### add your solver
