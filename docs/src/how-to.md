# How to...

[`WaterWaves1D.jl`](https://github.com/WaterWavesModels/WaterWaves1D.jl/) is meant to be versatile,
and integrating new blocks to the package is easy. If you ever do so, please do not hesitate to contact the
[developers](home.md#Developers) to either get help, or report on your advances.

## build your model

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

## build your initial data

The simplest way to build an initial data is to use the function [`Init`](@ref WaterWaves1D.Init), which takes as argument either
- a function `η` and a function `v` (in this order);
- an array of collocation points and two vectors representing `η(x)` and `v(x)` (in this order);
- a `mesh` (generated with [`Mesh`](@ref WaterWaves1D.Mesh)) and two vectors representing `η(mesh.x)` and `v(mesh.x)` (in this order).

Some relevant initial data (e.g. travelling waves) are built-in; see [the index section](index.md#Initial-data). You can build your own in the following lines.
```julia
struct Heap <: InitialData
	η
	v
	label :: String
	info  :: String

	function Heap(L)
		η=x->exp.(-(L*x).^2)
    v=x->zero(x)
		init = Init(η,v)
		label = "Heap"
		info = "Heap of water, with length L=$L."

		new( init.η,init.v,label,info  )
	end
end
```

The corresponding initial data can now be built as follows
```julia
using WaterWaves1D
# include your file
init = Heap(1)
```

## build your time solver

As an example, let us review how the explicit Euler solver, [`Euler`](@ref WaterWaves1D.Euler), is built.


In a dedicated file we write
```julia
struct Euler <: TimeSolver
    U1 :: Array
    label :: String

    function Euler( U :: Array; realdata=nothing )
        U1 = copy(U)
        if realdata==true
            U1 = real.(U1);
        end
        if realdata==false
            U1 = complex.(U1);
        end
        new( U1, "Euler" )
    end
end
```

Here, `Euler.U1` is a pre-allocated vector which can be used to speed-up calculations, and `Euler-label` is the string `"Euler"`, used for future references. The optional keyword argument `realdata` allows to specify the type of data which the solver will take as arguments: either complex or real vectors. 

We shall now add one method to the function `step!`, performing the explicit Euler step: that is replacing a vector `U` with `U+dt*f(U)` where `f` is provided by the model at stake.
```julia
export step!
function step!(solver :: Euler,
                model :: AbstractModel,
                U  ,
                dt )

    solver.U1 .= U        # allocate U to U1
    model.f!( solver.U1 ) # model.f!(U) replaces its argument U with f(U)
    U .+= dt .* solver.U1 # update U
end
```

## access to and manage your data

Once an initial-value problem `problem` has been solved (i.e. [numerically integrated](problems.md)), the raw data is stored in `problem.data.U`, which is an array whose elements correspond (in chronological order) to values at the different computed times, `problem.times.ts`. Physical data (say at final computed time) can be reconstructed using the model as follows: 
```julia
η,v,x = problem.model.mapfro(last(problem.data.U))
```
where  `η` and `v` are respectively the values of the surface deformation and of the derivative of the trace of the velocity potential at the surface, at collocation points `x`.

This procedure is carried out by the function [`solution`](@ref WaterWaves1D.solution), which allows in addition to perform some interpolations (making use of the otherwise helpful function [`interpolate`](@ref WaterWaves1D.interpolate)).

Using `(η,v,x)` one can compute other quantities such as the [mass, momentum, energy](background.md#Mass,-momentum,-energy); for instance for the purpose of testing how well these quantities are numerically perserved (when the quantities are first integrals of the considered model). Built-in functions [`mass`](@ref WaterWaves1D.mass), [`momentum`](@ref WaterWaves1D.momentum), [`energy`](@ref WaterWaves1D.energy) (and [`massdiff`](@ref WaterWaves1D.massdiff), [`momentumdiff`](@ref WaterWaves1D.momentumdiff), [`energydiff`](@ref WaterWaves1D.energydiff)) compute such quantities (or their variation).

## plot your data

Once `(η,v,x)` is obtained as above, producing plots is as simple as
```julia
using Plots
plot(x, [ η v ])
```
One can also plot the amplitude of [discrete Fourier coefficients](background.md#Pseudospectral-methods) (in semi-log scale) as follows:
```julia
using FFTW
k=fftshift(Mesh(x).k);fftη=fftshift(fft(η));fftv=fftshift(fft(v));
indices = (fftη .!=0) .& (fftv .!=0 )
plot(k[indices], [ abs.(fftη)[indices] abs.(fftv)[indices] ], yscale=:log10)
```

The convenient built-in function [`plot_solution`](@ref WaterWaves1D.plot_solution) produces such plots, allowing in addition interpolations or compression.

A similar function is [`plot_difference`](@ref WaterWaves1D.plot_difference) which allows to plot the difference between the outcome of several problems.

Finally, [`create_animation`](@ref WaterWaves1D.create_animation) allows to produce an animation with, for instance,
```julia
anim = create_animation(problem)
gif(anim, "my_animation.gif", fps=15)
```

## save and load your data
