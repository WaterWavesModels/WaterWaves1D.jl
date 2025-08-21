export KG

"""
    KG(param;kwargs)

Define an object of type `AbstractModel` in view of solving the initial-value problem for
the Klein-Gordon system 

# Argument
`param` is of type `NamedTuple` and must contain
- relaxation parameter `ϵ` (nonlinearity) and `δ` (dispersion);
- numerical parameters to construct the mesh of collocation points, if `mesh` is not provided as a keyword argument.

## Optional keyword arguments
- `mesh`: the mesh of collocation points. By default, `mesh = Mesh(param)`;
- `ktol`: tolerance of the Krasny filter (default is `0`, i.e. no filtering);
- `dealias`: dealiasing with Orszag rule `1-dealias/(dealias+2)` (default is `0`, i.e. no dealiasing);
- `label`: a label for future references;

# Return values
Generate necessary ingredients for solving an initial-value problem via `solve!`:
1. a function `KG.f!` to be called in explicit time-integration solvers;
2. a function `KG.mapto` which from `(η,v)` of type `InitialData` provides the raw data matrix on which computations are to be executed;
3. a function `KG.mapfro` which from such data matrix returns the Tuple of real vectors `(η,v,x)`, where
	- `η` is the values of surface deformation at collocation points `x`;
	- `v` is the derivative of the trace of the velocity potential at `x`;
4. additionally, a handy function `KG.mapfrofull` which from data matrix returns the Tuple of real vectors `(η,v,u,p,w)`, where
	- `u` corresponds to the layer-averaged horizontal velocity.
	- `p` corresponds to the relaxed (artificial) layer-averaged non-hydrostatic pressure;
	- `w` corresponds to the relaxed (artificial) layer-averaged vertical velocity.

"""
mutable struct KG <: AbstractModel

	label   :: String
	f!		:: Function
	mapto	:: Function
	mapfro	:: Function
	info    :: String

    function KG(param::NamedTuple; 
								mesh = Mesh(param),
								dealias	= 0,
								ktol	= 0,
								label	= nothing,
								c₀ = 0
								)
		# Set up
		δ 	= param.δ
		ϵ 	= param.ϵ

		if isnothing(label)
				label = "Klein-Gordon system"
		end


		# Print information
		info = "$label model.\n"
		info *= "├─Relaxation parameter ϵ=$ϵ.\n"
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
		x 	= mesh.x
		∂ₓ	=  1im * k

		F 	= sqrt.(1 .+ abs.(δ * k).^2)


		if dealias == 0
			Π⅔ 	= ones(size(k)) # no dealiasing (Π⅔=Id)
		else
			K = (mesh.kmax-mesh.kmin)/(2+dealias)
			Π⅔ 	= abs.(k) .<= K # Dealiasing low-pass filter
		end
		u = zeros(Complex{Float64}, mesh.N)
		fftv, fftw = (similar(u),).*ones(2)

		# Evolution equations are ∂t U = f(U)
		function f!(U;ϵ=ϵ)
				fftv .= U[1]; fftw .= U[2] ;
	
				U[1] .= -1/ϵ * F.* fftw
				U[2] .=  1/ϵ * F.* fftv .- c₀*(δ*∂ₓ./F).^2 .* (∂ₓ.*fftw)

				for u in U u[ abs.(u).< ktol ].=0 end
		end

		# Build raw data from physical data.
		# Discrete Fourier transform with, possibly, dealiasing and Krasny filter.
		function mapto(data::InitialData)
			
			U = [Π⅔.*fft(data.η(x)), Π⅔.*fft(data.v(x))]

			for u in U u[ abs.(u).< ktol ].=0 end
			return U
		end

		# Reconstruct variables from raw data
		# Return `(η,v,x)`, where
		# - `u` is the surface deformation;
		# - `v` is the derivative of the trace of the velocity potential;
		# - `x` is the vector of collocation points
		function mapfro(U)
			real(ifft(U[1])),real(ifft(U[2])),mesh.x
		end

		new(label, f!, mapto, mapfro, info )
    end


end
