export Burgers

"""
    Burgers(param;kwargs)

Define an object of type `AbstractModel` which is a Burgers model for three-scale singular problems

# Argument
`param` is of type `NamedTuple` and must contain
- the nonlinearity parameter `ϵ`;
- numerical parameters to construct the mesh of collocation points, if `mesh` is not provided as a keyword argument.

## Optional keyword arguments
- `mesh`: the mesh of collocation points. By default, `mesh = Mesh(param)`;
- `ktol`: tolerance of the Krasny filter (default is `0`, i.e. no filtering);
- `dealias`: dealiasing with Orlicz rule `1-dealias/(dealias+2)` (default is `0`, i.e. no dealiasing);
- `label`: a label for future references (default is `"Burgers");

# Return values
Generate necessary ingredients for solving an initial-value problem via `solve!`:
1. a function `Burgers.f!` to be called in explicit time-integration solvers;
2. a function `Burgers.mapto` which from `(η,v)` of type `InitialData` provides the raw data matrix on which computations are to be executed;
3. a function `Burgers.mapfro` which from such data matrix returns the Tuple of real vectors `(h,u,x)`, where
	- `h` is the values of the static component at `x`;
	- `u` is the values of the real part of solutions at collocation points `x`;
4. a function `Burgers.mapfrofull` which from such data matrix returns the Tuple of real vectors `(h,u,x)`, where
	- `h` is the values of the static component at `x`;
	- `u` is the values of solutions at collocation points `x`;

"""
mutable struct Burgers <: AbstractModel

	label   :: String
	f!		:: Function
	mapto	:: Function
	mapfro	:: Function
	info    :: String

    function Burgers(param::NamedTuple; 
								mesh = Mesh(param),
								str	= true,
								dealias = 0,
								ktol	= 0,
								label	= nothing
								)
		# Set up
		ϵ	= param.ϵ

		if isnothing(label)
				label = "Burgers"
		end


		# Print information
		info = "$label model.\n"
		info *= "├─nonlinarity parameter ϵ=$ϵ.\n"

		if dealias == 0
			info *= "├─No dealiasing. "
		else
			info *= "├─Dealiasing with Orszag's rule adapted to power $(dealias + 1) nonlinearity. "
		end
		if ktol == 0
			info *= "No Krasny filter. "
		else
			info *= "Krasny filter with tolerance $ktol."
		end
		info *= "\nDiscretized with $(mesh.N) collocation points on [$(mesh.xmin), $(mesh.xmax)]."

		# Pre-allocate useful data
		k = mesh.k
		x = mesh.x

		∂ₓ	=  1im * k

		if dealias == 0
			Π⅔ 	= ones(size(k)) # no dealiasing (Π⅔=Id)
		else
			K = (mesh.kmax-mesh.kmin)/3
            Π⅔      =min.(( abs.(k) .<= K).*(-abs.(k)/dealias.+K/dealias),1)		
        end
		u = zeros(Complex{Float64}, mesh.N)


		# Evolution equations are ∂t U = f(U)
			function f!(U)
				U .= -ϵ/2 * Π⅔.*∂ₓ.*fft( ifft( U ).^2 )
				U[ abs.(U).< ktol].=0
			end

		# Build raw data from physical data.
		# Discrete Fourier transform with, possibly, dealiasing and Krasny filter.
		function mapto(data::InitialData)
			u .= fft(data.v(x)) 
			u[abs.(u).< ktol ].=0
			U = reshape(Π⅔.*u,length(u),1)
			return U
		end

		# Reconstruct physical variables from raw data
		# Return `(u,x)`, where
		# - `u` is the values of the real part of solutions at collocation points `x`;
		# - `x` is the vector of collocation points
		function mapfro(U)
			u .= ifft(U)

			return real(u),real(u),mesh.x 
		end

        new(label, f!, mapto, mapfro, info )
    end


end
