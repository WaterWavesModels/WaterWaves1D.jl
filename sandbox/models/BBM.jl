export BBM

"""
    BBM(param;kwargs)

Define an object of type `AbstractModel` in view of solving the initial-value problem for
the BBM equation

# Argument
`param` is of type `NamedTuple` and must contain
- dimensionless parameters `ϵ` (nonlinearity) and `μ` (dispersion);
- numerical parameters to construct the mesh of collocation points, if `mesh` is not provided as a keyword argument.

## Optional keyword arguments
- `KdV`: if `true` (default is `false`), compute the standard KdV equations instead (see `KdV(param;kwargs)`);
- `mesh`: the mesh of collocation points. By default, `mesh = Mesh(param)`;
- `ktol`: tolerance of the low-pass Krasny filter (default is `0`, i.e. no filtering);
- `dealias`: dealiasing with Orszag rule `1-dealias/(dealias+2)` (default is `0`, i.e. no dealiasing);
- `label`: a label for future references (default is `"BBM-Boussinesq"`);


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
								mesh = Mesh(param),
								dealias = 0,
								ktol	= 0,
								label 	= nothing
								)

		# Set up
		δ 	= param.δ
		if isnothing(label) label = "BBM" end
		
		# Print information
		info = "$label model.\n"
		info *= "├─Shallowness parameter δ=$δ.\n"
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
		info *= "\nDiscretized with $(mesh.N) collocation points on [$(mesh.xmin), $(mesh.xmax)]."

		# Pre-allocate data
		x 	= mesh.x
		k = mesh.k
		∂ₓ	=  1im * k
		F 	= ∂ₓ./(1 .+ abs.(δ * k).^2)
		
		if dealias == 0
			Π⅔ 	= ones(size(k)) # no dealiasing (Π⅔=Id)
		else
			K = (mesh.kmax-mesh.kmin)/(2+dealias)
			Π⅔ 	= abs.(k) .<= K # Dealiasing low-pass filter
		end
		u = zeros(Float64, mesh.N)
		fftu = zeros(Complex{Float64}, mesh.N)

		# Evolution equations are ∂t U = f(U)
		function f!(U)

		    fftu .= U;u.=real.(ifft(fftu))

		   	U .= -F.*(fftu+1/2*Π⅔.*fft(u.^2))
			
			U[ abs.(u).< ktol ].=0

		end

		# Build raw data from physical data.
		# Discrete Fourier transform of the suitable variables 
		# with, possibly, dealiasing and Krasny filter.
		function mapto(data::InitialData)
			U = fft(data.v(x)) 
			U[ abs.(U).< ktol ].=0 
			return U
		end

		# Reconstruct physical variables from raw data
		# Return `(u,x)`, where
		# - `u` is the solution;
		# - `x` is the vector of collocation points
		function mapfro(U)
			real(ifft(U)),real(ifft(U)),mesh.x
		end

        new(label, f!, mapto, mapfro, info )
    end
end
