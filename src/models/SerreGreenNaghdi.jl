export SerreGreenNaghdi
using FFTW,LinearAlgebra,LinearMaps,IterativeSolvers

"""
    SerreGreenNaghdi(param;kwargs)

Define an object of type `AbstractModel` in view of solving the initial-value problem for
the Serre-Green-Naghdi model

# Argument
`param` is of type `NamedTuple` and must contain
- dimensionless parameters `ϵ` (nonlinearity) and `μ` (dispersion);
- numerical parameters to construct the mesh of collocation points as `mesh = Mesh(param)`.

## Keywords
- `iterative`: solve the elliptic problem through GMRES if `true`, LU decomposition if `false` (default is `true`);
- `precond`: use a (left) preconditioner for GMRES if `true` (default), choose `precond` as the preconditioner if provided;
- `gtol`: relative tolerance of the GMRES algorithm (default is `1e-14`);
- `restart`: the corresponding option of the GMRES algorithm (default is `100`);
- `maxiter`: the corresponding option of GMRES (default is `nothing`);
- `ktol`: tolerance of the Krasny filter (default is `0`, i.e. no filtering);
- `dealias`: dealiasing with Orlicz rule `1-dealias/(dealias+2)` (default is `0`, i.e. no dealiasing);
- `verbose`: prints information if `true` (default is `true`).

# Return values
Generate necessary ingredients for solving an initial-value problem via `solve!` and in particular
1. a function `SerreGreenNaghdi.f!` to be called in the time-integration solver;
2. a function `SerreGreenNaghdi.mapto` which from `(η,v)` of type `InitialData` provides the raw data matrix on which computations are to be executed;
3. a function `SerreGreenNaghdi.mapfro` which from such data matrix returns the Tuple of real vectors `(η,v)`, where

    - `η` is the surface deformation;
    - `v` is the derivative of the trace of the velocity potential;
4. additionally, a handy function `SerreGreenNaghdi.mapfrofull` which from data matrix returns the Tuple of real vectors `(η,v,u)`, where

	- `u` corresponds to the layer-averaged velocity.

"""
mutable struct SerreGreenNaghdi <: AbstractModel

	label   :: String
	f!		:: Function
	mapto	:: Function
	mapfro	:: Function
	mapfrofull	:: Function
	param	:: NamedTuple
	kwargs  :: NamedTuple

    function SerreGreenNaghdi(param::NamedTuple;dealias=0,ktol=0,iterate=true,gtol=1e-14,precond=true,restart=nothing,maxiter=nothing,verbose=true)
		m=WhithamGreenNaghdi(param;SGN=true,dealias=dealias,ktol=ktol,iterate=iterate,gtol=gtol,precond=precond,restart=restart,maxiter=maxiter,verbose=verbose)
        new(m.label, m.f!, m.mapto, m.mapfro, m.mapfrofull, m.param, m.kwargs)
    end
end
