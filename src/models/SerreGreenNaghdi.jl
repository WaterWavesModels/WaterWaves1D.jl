export SerreGreenNaghdi

"""
    SerreGreenNaghdi(param; kwargs...)

Define an object of type `AbstractModel` in view of solving the initial-value problem for
the Serre-Green-Naghdi model ([Serre](@cite Serre1953), [Su and Gardner](@cite Su1969), [Green and Naghdi](@cite Green1976)).

# Argument
`param` is of type `NamedTuple` and must contain
- dimensionless parameters `ϵ` (nonlinearity) and `μ` (dispersion);
- numerical parameters to construct the mesh of collocation points, if `mesh` is not provided as a keyword argument.

## Optional keyword arguments
- `mesh`: the mesh of collocation points. By default, `mesh = Mesh(param)`;
- `iterative`: solve the elliptic problem through GMRES if `true`, LU decomposition if `false` (default is `true`);
- `precond`: use a (left) preconditioner for GMRES if `true` (default), choose `precond` as the preconditioner if provided;
- `gtol`: relative tolerance of the GMRES algorithm (default is `1e-14`);
- `restart`: the corresponding option of the GMRES algorithm (default is `100`);
- `maxiter`: the corresponding option of GMRES (default is `nothing`);
- `ktol`: tolerance of the Krasny filter (default is `0`, i.e. no filtering);
- `dealias`: dealiasing with Orszag rule `1-dealias/(dealias+2)` (default is `0`, i.e. no dealiasing);
- `label`: a label for future references (default is `"Green-Naghdi"`);

# Return values
Generate necessary ingredients for solving an initial-value problem via `solve!`:
1. a function `SerreGreenNaghdi.f!` to be called in explicit time-integration solvers;
2. a function `SerreGreenNaghdi.mapto` which from `(η,v)` of type `InitialData` provides the raw data matrix on which computations are to be executed;
3. a function `SerreGreenNaghdi.mapfro` which from such data matrix returns the Tuple of real vectors `(η,v,x)`, where
    - `η` is the values of surface deformation at collocation points `x`;
    - `v` is the derivative of the trace of the velocity potential at `x`;
4. additionally, a handy function `SerreGreenNaghdi.mapfrofull` which from data matrix returns the Tuple of real vectors `(η,v,u)`, where
    - `u` corresponds to the layer-averaged velocity.

"""
mutable struct SerreGreenNaghdi <: AbstractModel

	label   :: String
	f!		:: Function
	mapto	:: Function
	mapfro	:: Function
	mapfrofull	:: Function
	info	:: String

    function SerreGreenNaghdi(param::NamedTuple;
							mesh = Mesh(param),
							dealias = 0,
							ktol	= 0,
							iterate	= true,
							gtol	= 1e-14,
							precond	= true,
							restart	= nothing,
							maxiter	= nothing,
							label	= "Green-Naghdi"
							)

		m=WhithamGreenNaghdi(param;SGN=true,
							mesh=mesh,
							dealias=dealias,
							ktol=ktol,
							iterate=iterate,
							gtol=gtol,
							precond=precond,
							restart=restart,
							maxiter=maxiter,
							label=label
							)

		new(m.label, m.f!, m.mapto, m.mapfro, m.mapfrofull, m.info )
    end
end
