export Toy

"""
    Toy(param;kwargs)

Define an object of type `AbstractModel` which is a toy model for three-scale singular problems

# Argument
`param` is of type `NamedTuple` and must contain
- the relaxation parameters `a` and `δ`;
- numerical parameters to construct the mesh of collocation points, if `mesh` is not provided as a keyword argument.

## Optional keyword arguments
- `mesh`: the mesh of collocation points. By default, `mesh = Mesh(param)`;
- `h`: the static component function. By default, `h = sin`;
- `ktol`: tolerance of the Krasny filter (default is `0`, i.e. no filtering);
- `dealias`: dealiasing with Orszag rule `1-dealias/(dealias+2)` (default is `0`, i.e. no dealiasing);
- `label`: a label for future references (default is `"toy");

# Return values
Generate necessary ingredients for solving an initial-value problem via `solve!`:
1. a function `toy.f!` to be called in explicit time-integration solvers;
2. a function `toy.mapto` which from `(η,v)` of type `InitialData` provides the raw data matrix on which computations are to be executed;
3. a function `toy.mapfro` which from such data matrix returns the Tuple of real vectors `(h,u,x)`, where
	- `h` is the values of the static component at `x`;
	- `u` is the values of the real part of solutions at collocation points `x`;
4. a function `toy.mapfrofull` which from such data matrix returns the Tuple of real vectors `(h,u,x)`, where
	- `h` is the values of the static component at `x`;
	- `u` is the values of solutions at collocation points `x`;

"""
mutable struct Toy <: AbstractModel

	label   :: String
	f!		:: Function
	mapto	:: Function
	mapfro	:: Function
	mapfrofull	:: Function
	info    :: String

    function Toy(param::NamedTuple; 
								mesh = Mesh(param),
								str	= true,
								dealias = 0,
								ktol	= 0,
								label	= nothing
								)
		# Set up
		a	= param.a
		δ 	= param.δ

		if isnothing(label)
				label = "toy"
		end


		# Print information
		info = "$label model.\n"
		info *= "├─Relaxation parameter a=$a.\n"
		info *= "├─Shallowness parameter δ=$δ.\n"

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
		Jδ = δ * ∂ₓ .+ 1im

		if dealias == 0
			Π⅔ 	= ones(size(k)) # no dealiasing (Π⅔=Id)
		else
			K = (mesh.kmax-mesh.kmin)/(2+dealias)
			Π⅔ 	= abs.(k) .<= K # Dealiasing low-pass filter
		end
		h = zeros(Complex{Float64}, mesh.N)
		u = zeros(Complex{Float64}, mesh.N)


		# Evolution equations are ∂t U = f(U)
			function f!(U)
				if str == true
					U[:,2] .= -a * Π⅔.*fft( ifft(Jδ  .* U[:,2] ) ./ U[:,1] )
				else
					U[:,2] .= -a * Π⅔.*fft( ifft(Jδ  .* U[:,2] ) ./ ( 1 .+ abs.(ifft(U[:,2])).^2 ) )
				end
				U[ abs.(U[:,2]).< ktol, 2 ].=0
				U[:,1] .= zero(U[:,1])
			end
			#  function f!(U)
			# 	U[:,2] .= -a * Π⅔.*fft( ifft(Jδ  .* U[:,2] ) .* ( 1 .+ ifft(U[:,2]).^2 ) )
			# 	U[ abs.(U[:,2]).< ktol, 2 ].=0
			# 	U[:,1] = zero(U[:,1])
			# end

		# Build raw data from physical data.
		# Discrete Fourier transform with, possibly, dealiasing and Krasny filter.
		function mapto(data::InitialData)
			h .= 1 .+ data.η(x)
			u .= fft(data.v(x)) 
			u[abs.(u).< ktol ].=0
			
			U = [h Π⅔ .* u]
			return U
		end

		# Reconstruct physical variables from raw data
		# Return `(h,u,x)`, where
		# - `h` is the values of the static component at `x`;
		# - `u` is the values of the real part of solutions at collocation points `x`;
		# - `x` is the vector of collocation points
		function mapfro(U)
			h .= U[:,1]
			u .= ifft(U[:,2])

			return real(h),real(u),mesh.x 
		end
		
		# Reconstruct physical variables from raw data
		# Return `(h,u,x)`, where
		# - `h` is the values of the static component at `x`;
		# - `u` is the values of the real part of solutions at collocation points `x`;
		# - `x` is the vector of collocation points
		function mapfrofull(U)
			h .= U[:,1]
			u .= ifft(U[:,2])

			return real(h),u,mesh.x 
		end

        new(label, f!, mapto, mapfro, mapfrofull, info )
    end


end
