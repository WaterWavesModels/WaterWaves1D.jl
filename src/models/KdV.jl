export KdV

@doc raw"""
    KdV(param; kwargs...)

Define an object of type `AbstractModel` in view of solving the initial-value problem for
two uncoupled [KortewegDeVries1895](@citet) equations:
```math
∂_t u_\pm   \mp ∂_x u_\pm \mp \frac{\mu}{6} ∂_x^3 u_\pm \mp \frac{3ϵ}{4} ∂_x (u_\pm^2) = 0 ,
```
where  ``u_\pm`` are right- and left-going waves, combined to produce the physical solutions.

See [Lannes2013](@citet) and [Emerald2021b](@citet) for formulas and a full justification.

# Argument

`param` is of type `NamedTuple` and must contain
- dimensionless parameters `ϵ` (nonlinearity) and `μ` (dispersion);
- numerical parameters to construct the mesh of collocation points, if `mesh` is not provided as a keyword argument.

## Optional keyword arguments

- `improved_initial_data`: if `true` (default), improves the naive (first-order) decomposition into right-going and left-going wave;
- `mesh`: the mesh of collocation points. By default, `mesh = Mesh(param)`;
- `ktol`: tolerance of the low-pass Krasny filter (default is `0`, i.e. no filtering);
- `dealias`: dealiasing with Orszag rule `1-dealias/(dealias+2)` (default is `0`, i.e. no dealiasing);
- `label`: a label for future references (default is `"KdV"`).


# Return values
Generate necessary ingredients for solving an initial-value problem via `solve!`:
1. a function `KdV.f!` to be called in explicit time-integration solvers;
Additionnally, a vector `Whitham.D` and a function `Whitham.f!` for exponential solvers;
2. a function `KdV.mapto` which from `(η,v)` of type `InitialData` provides the raw data matrix on which computations are to be executed.
3. a function `KdV.mapfro` which from such data matrix returns the Tuple of real vectors `(η,v,x)`, where
  - `η` is the values of surface deformation at collocation points `x`;
  - `v` is the derivative of the trace of the velocity potential at `x`.

"""
mutable struct KdV <: AbstractModel

	label   :: String
	f!		:: Function
	D		:: AbstractArray
	g!		:: Function
	mapto	:: Function
	mapfro	:: Function
	info 	:: String

    function KdV(param::NamedTuple;
								improved_initial_data=true,
								mesh = Mesh(param),
								dealias = 0,
								ktol	= 0,
								label 	= nothing
								)

		m=Whitham(param;KdV=true,
							improved_initial_data=improved_initial_data,
							mesh=mesh,
							dealias=dealias,ktol=ktol,
							label=label)

		new(m.label, m.f!, m.D, m.g!, m.mapto, m.mapfro, m.info)
    end
end
