export Boussinesq

@doc raw"""
    Boussinesq(param; kwargs...)

Define an object of type `AbstractModel` in view of solving the initial-value problem for
`abcd`-Boussinesq models proposed in [BonaChenSaut2002](@citet):
```math
  \left\{\begin{array}{l}
  ∂_tη+∂_x((1+ϵη)u)+a∂_x^3u-b∂_x^2∂_tη=0,\\[1ex]
  ∂_tu+∂_xη+ϵu∂_xu+c∂_x^3η-d∂_x^2∂_tu=0,
  \end{array}\right.
```
where ``a,b,c,d`` can be freely chosen as long as they satisfy ``a+b+c+d=1/3``.

The function is restricted to the so-called Hamiltonian case: ``d=b``, and ``c=0``.
Indeed, in that case, one may genuinely interpret ``v=u-d∂_x^2u`` as an approximation to ``∂_xψ``, the derivative of the trace of the velocity potential at the surface, and write the system as evolution equations for the variables ``(η,v)``:
```math
  \left\{\begin{array}{l}
  ∂_tη+∂_x((1+a∂_x^2)(1-b∂_x^2)^{-2}v + ϵ(1-b∂_x^2)^{-1} (η (1-b∂_x^2)^{-1}v))=0,\\[1ex]
  ∂_tv+∂_xη+\tfrac{ϵ}{2}∂_x(((1-b∂_x^2)^{-1}v)^2) =0.
  \end{array}\right.
```

# Argument
`param` is of type `NamedTuple` and must contain
- dimensionless parameters `ϵ` (nonlinearity) and `μ` (dispersion);
- numerical parameters to construct the mesh of collocation points, if `mesh` is not provided as a keyword argument.

## Optional keyword arguments
- two parameters `a` (default is `-1/3`) and `b` (default is `+1/3`) which determine the model solved. You need `a+2*b=1/3` for validity as a long wave model (without surface tension).
- `mesh`: the mesh of collocation points. By default, `mesh = Mesh(param)`;
- `ktol`: tolerance of the low-pass Krasny filter (default is `0`, i.e. no filtering);
- `dealias`: dealiasing with Orszag rule `1-dealias/(dealias+2)` (default is `0`, i.e. no dealiasing);
- `label`: a label for future references (default is `"Boussinesq"`);

# Return values
Generate necessary ingredients for solving an initial-value problem via `solve!`:
1. a function `Boussinesq.f!` to be called in explicit time-integration solvers;
2. a function `Boussinesq.mapto` which from `(η,v)` of type `InitialData` provides the raw data matrix on which computations are to be executed;
3. a function `Boussinesq.mapfro` which from such data matrix returns the Tuple of real vectors `(η,v,x)`, where
    - `η` is the values of surface deformation at collocation points `x`;
    - `v` is the derivative of the trace of the velocity potential at `x`.

"""
mutable struct Boussinesq <: AbstractModel

	label 	:: String
	f!		:: Function
	mapto	:: Function
	mapfro	:: Function
	info	:: String

    function Boussinesq(param::NamedTuple;
						mesh = Mesh(param),
						a=-1//3,b=1//3,
						dealias=0,ktol=0,
						label="Boussinesq"
						)

		m=WhithamBoussinesq(param;Boussinesq=true,
							mesh=mesh,
							a=a,b=b,
							dealias=dealias,ktol=ktol,
							label=label)

		new(m.label, m.f!, m.mapto, m.mapfro, m.info)
    end
end
