export SaintVenant

"""
    SaintVenant(param;kwargs)

Define an object of type `AbstractModel` in view of solving the initial-value problem for
Saint-Venant (or shallow water) model.

# Argument
`param` is of type `NamedTuple` and must contain
- the dimensionless parameter `ϵ` (nonlinearity);
- numerical parameters to construct the mesh of collocation points as `mesh = Mesh(param)`.

## Optional keyword arguments
- `ktol`: tolerance of the low-pass Krasny filter (default is `0`, i.e. no filtering);
- `dealias`: dealiasing with Orlicz rule `1-dealias/(dealias+2)` (default is `0`, i.e. no dealiasing);
- `label`: a label for future references (default is `"Saint-Venant"`);

# Return values
Generate necessary ingredients for solving an initial-value problem via `solve!`:
1. a function `SaintVenant.f!` to be called in explicit time-integration solvers;
2. a function `SaintVenant.mapto` which from `(η,v)` of type `InitialData` provides the raw data matrix on which computations are to be executed;
3. a function `SaintVenant.mapfro` which from such data matrix returns the Tuple of real vectors `(η,v)`, where
    - `η` is the surface deformation;
    - `v` is the derivative of the trace of the velocity potential.

"""
mutable struct SaintVenant <: AbstractModel

	label 	:: String
	f!		:: Function
	mapto	:: Function
	mapfro	:: Function
	info	:: String

    function SaintVenant(param::NamedTuple;
						dealias=0,ktol=0,
						label="Saint-Venant"
						)

		m=WhithamBoussinesq(merge(param,(μ=1,));Boussinesq=true,
							a=0,b=0,
							dealias=dealias,ktol=ktol,
							label=label)

		info = "Saint-Venant model.\n"
		info *= "├─Nonlinearity parameter ϵ=$(param.ϵ).\n"
		if dealias == 0
			info *= "└─No dealiasing. "
		else
			info *= "└─Dealiasing with Orszag's rule adapted to power $(dealias + 1) nonlinearity. "
		end
		if ktol == 0
			info *= "No Krasny filter. "
		else
			info *= "Krasny filter with tolerance $ktol."
		end
		mesh=Mesh(param)
		info *= "\nDiscretized with $(mesh.N) collocation points on [$(mesh.xmin), $(mesh.xmax)]."


		new(m.label, m.f!, m.mapto, m.mapfro, info)
    end
end
