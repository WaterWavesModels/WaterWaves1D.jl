export BBM

"""
    BBM(param; kwargs...)

Define an object of type `AbstractModel` in view of solving the initial-value problem for
two uncoupled Benjamin-Bona-Mahony equations.

# Argument
`param` is of type `NamedTuple` and must contain
- dimensionless parameters `ϵ` (nonlinearity) and `μ` (dispersion);
- numerical parameters to construct the mesh of collocation points, if `mesh` is not provided as a keyword argument.

## Optional keyword arguments
- `improved_initial_data`: if `true` (default), improves the naive (first-order) decomposition into right-going and left-going wave;
- `mesh`: the mesh of collocation points. By default, `mesh = Mesh(param)`;
- `ktol`: tolerance of the low-pass Krasny filter (default is `0`, i.e. no filtering);
- `dealias`: dealiasing with Orszag rule `1-dealias/(dealias+2)` (default is `0`, i.e. no dealiasing);
- `label`: a label for future references (default is `"BBM"`).


# Return values
Generate necessary ingredients for solving an initial-value problem via `solve!`:
1. a function `BBM.f!` to be called in explicit time-integration solvers;
2. a function `BBM.mapto` which from `(η,v)` of type `InitialData` provides the raw data matrix on which computations are to be executed.
3. a function `BBM.mapfro` which from such data matrix returns the Tuple of real vectors `(η,v,x)`, where
  - `η` is the values of surface deformation at collocation points `x`;
  - `v` is the derivative of the trace of the velocity potential at `x`.

"""
mutable struct BBM <: AbstractModel

	label   :: String
	f!		:: Function
	mapto	:: Function
	mapfro	:: Function
	info 	:: String

    function BBM(param::NamedTuple;
								improved_initial_data=true,
								mesh = Mesh(param),
								dealias = 0,
								ktol	= 0,
								label 	= nothing
								)

		m=Whitham(param;BBM=true,
							improved_initial_data=improved_initial_data,
							mesh=mesh,
							dealias=dealias,ktol=ktol,
							label=label)

		new(m.label, m.f!, m.mapto, m.mapfro, m.info)
    end
end
