export SaintVenant2D

"""
    SaintVenant2D(param;kwargs)

Define an object of type `AbstractModel` in view of solving the initial-value problem for
Saint-Venant (or shallow water) model.

# Argument
`param` is of type `NamedTuple` and must contain
- the dimensionless parameter `ϵ` (nonlinearity);
- numerical parameters to construct the mesh of collocation points, if `mesh` is not provided as a keyword argument.

## Optional keyword arguments
- `mesh`: the mesh of collocation points. By default, `mesh = Mesh(param)`;
- `ktol`: tolerance of the low-pass Krasny filter (default is `0`, i.e. no filtering);
- `dealias`: no dealisasing if set to `0` or `false` (default), standard 3/2 Orszag rule if set to `1` or `true`, otherwise the value sets additionnally a maximal slope of the dealiasing symbol (`2/dealias` modes are affected);
- `label`: a label for future references (default is `"Saint-Venant"`);

# Return values
Generate necessary ingredients for solving an initial-value problem via `solve!`:
1. a function `SaintVenant2D.f!` to be called in explicit time-integration solvers;
2. a function `SaintVenant2D.mapto` which from `(η,v)` of type `InitialData` provides the raw data matrix on which computations are to be executed;
3. a function `SaintVenant2D.mapfro` which from such data matrix returns the Tuple of real vectors `(η,v,x)`, where
    - `η` is the values of surface deformation at collocation points `x`;
    - `v` is the derivative of the trace of the velocity potential at `x`.

"""
mutable struct SaintVenant2D <: AbstractModel

	label 	:: String
	f!		:: Function
	mapto	:: Function
	mapfro	:: Function
	info	:: String

    function SaintVenant2D(param::NamedTuple;
						mesh = Mesh(param),
						dealias=0,ktol=0,
						hamiltonian=false,
						label="Saint-Venant"
						)

		# Set up
		ϵ 	= param.ϵ

		if hamiltonian 
			info = "2D Saint-Venant model (hamiltonian).\n"
		else
			info = "2D Saint-Venant model.\n"
		end
		info *= "├─Nonlinearity parameter ϵ=$(param.ϵ).\n"
		if dealias == 0
			info *= "└─No dealiasing. "
		else
			info *= "├─Dealiasing with Orszag's rule adapted to power 2 nonlinearity: \n"
			if dealias == 1
				info *= "└─Sharp cut-off. "
			else
				info *= "└─Lipschitz cut-off. "
			end
		end
		if ktol == 0
			info *= "No Krasny filter. "
		else
			info *= "Krasny filter with tolerance $ktol."
		end
		info *= "\nDiscretized with $(mesh.N) collocation points on [$(mesh.xmin), $(mesh.xmax)]."

		# Pre-allocate useful data
		kx=mesh.k;ky=mesh.k;x=mesh.x;y=mesh.x;
		∂y = 1im * ky' .* one.( x)
		∂x = one.(y') .* 1im .* kx


		if dealias == 0
			Π⅔ 	= ones(mesh.N,mesh.N) # no dealiasing (Π⅔=Id)
			Π 	= ones(mesh.N,mesh.N) # no dealiasing (Π⅔=Id)
		else
			K = (mesh.kmax-mesh.kmin)/3
			Π⅔ 	= (abs.(mesh.k) .<= K)'.*(abs.(mesh.k) .<= K) # Dealiasing low-pass filter
			P   = min.(( abs.(mesh.k) .<= K).*( (-abs.(mesh.k) .+ K)*dealias/mesh.dk ),1)	
			#Π  =min.(( abs.(mesh.k) .<= K).*( (-abs.(mesh.k) .+ K)*dealias/mesh.dk .+1/2 ),1)	
			Π = P'.*P		
		end

		
		η = zeros(Float64, (mesh.N,mesh.N))
        vx = zeros(Float64, (mesh.N,mesh.N))
		vy = zeros(Float64, (mesh.N,mesh.N))
		fftη = zeros(Complex{Float64}, (mesh.N,mesh.N))
        fftvx = zeros(Complex{Float64}, (mesh.N,mesh.N))
		fftvy = zeros(Complex{Float64}, (mesh.N,mesh.N))

		# Evolution equations are ∂t U = f(U)
		function f!(U)
			fftη .= U[:,:,1]
			fftvx .= U[:,:,2]
			fftvy .= U[:,:,3]

			η .= real(ifft(U[:,:,1]))
			vx .= real(ifft(U[:,:,2]))
			vy .= real(ifft(U[:,:,3]))

			U[:,:,1] .= -∂x.*( fftvx .+ ϵ * Π.* fft( η .* vx) ) +
						-∂y.*( fftvy .+ ϵ * Π.* fft( η .* vy) )
			if hamiltonian
				U[:,:,2] .= -∂x.*( fftη .+ ϵ/2 * Π.* fft( vx.^2+vy.^2 ) )
				U[:,:,3] .= -∂y.*( fftη .+ ϵ/2 * Π.* fft( vx.^2+vy.^2 ) )
			else
				U[:,:,2] .= -∂x.*( fftη ) .- ϵ * Π.* fft( vx.* ifft(∂x.*fftvx) + vy.* ifft(∂y.*fftvx) )
				U[:,:,3] .= -∂y.*( fftη ) .- ϵ * Π.* fft( vx.* ifft(∂x.*fftvy) + vy.* ifft(∂y.*fftvy) )
			end		
			U[ abs.(U).< ktol ].=0
		end

		# Build raw data from physical data.
		# Discrete Fourier transform with, possibly, dealiasing and Krasny filter.
		function mapto(data::InitialData)
			U = cat(Π⅔.* fft(data.η(x,y)), Π⅔.*fft(data.vx(x,y)), Π⅔.*fft(data.vy(x,y)),dims=3)
			#U = [Π⅔.* fft(data.η(x,y)), Π⅔.*fft(data.vx(x,y)), Π⅔.*fft(data.vy(x,y))]

			U[abs.(U).<ktol].=0
			# U[1][ abs.(U[1]).< ktol ].=0
			# U[2][ abs.(U[2]).< ktol ].=0
			# U[3][ abs.(U[3]).< ktol ].=0
			return U
		end

		# Reconstruct physical variables from raw data
		# Return `(η,v,x)`, where
		# - `η` is the surface deformation;
		# - `v` is the derivative of the trace of the velocity potential;
		# - `x` is the vector of collocation points
		function mapfro(U)
			real(ifft(U[:,:,1])),real(ifft(U[:,:,2])),real(ifft(U[:,:,3])),x,y
		end

		new(label, f!, mapto, mapfro, info)
    end
end
