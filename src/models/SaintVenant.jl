export SaintVenant,SaintVenant_fast

"""
    SaintVenant(param;kwargs)

Define an object of type `AbstractModel` in view of solving the initial-value problem for
Saint-Venant (or shallow water) model.

# Argument
`param` is of type `NamedTuple` and must contain
- the dimensionless parameter `ϵ` (nonlinearity);
- numerical parameters to construct the mesh of collocation points, if `mesh` is not provided as a keyword argument.

## Optional keyword arguments
- `mesh`: the mesh of collocation points. By default, `mesh = Mesh(param)`;
- `ktol`: tolerance of the low-pass Krasny filter (default is `0`, i.e. no filtering);
- `dealias`: no dealisasing if set to `0` or `false` (default), otherwise `1/(3*dealias)` modes are set to `0` (corresponding to standard 2/3 Orszag rule if `dealias` is set to `1` or `true`);
- `smooth`: A smooth low-pass filter (whose scaling is defined by ) if set to `0` or `false` (default), otherwise only `2/(3*dealias)*(1-smooth/2)` modes are kept untouched;
- `label`: a label for future references (default is `"Saint-Venant"`);

# Return values
Generate necessary ingredients for solving an initial-value problem via `solve!`:
1. a function `SaintVenant.f!` to be called in explicit time-integration solvers;
2. a function `SaintVenant.mapto` which from `(η,v)` of type `InitialData` provides the raw data matrix on which computations are to be executed;
3. a function `SaintVenant.mapfro` which from such data matrix returns the Tuple of real vectors `(η,v,x)`, where
    - `η` is the values of surface deformation at collocation points `x`;
    - `v` is the layer-averaged velocity (or the derivative of the trace of the velocity potential) at `x`.

"""
mutable struct SaintVenant <: AbstractModel

	label 	:: String
	f!		:: Function
	mapto	:: Function
	mapfro	:: Function
	info	:: String

    function SaintVenant(param::NamedTuple;
						mesh = Mesh(param),
						dealias=0,smooth=false,
						ktol=0,
						label="Saint-Venant"
						)

		# Set up
		ϵ 	= param.ϵ

		info = "Saint-Venant model.\n"
		info *= "├─Nonlinearity parameter ϵ=$(param.ϵ).\n"
		if dealias == 0
			info *= "└─No dealiasing. "
		elseif dealias == 1
			info *= "├─Dealiasing with Orszag's rule adapted to power 2 nonlinearity. \n"
		else
			info *= "├─Dealiasing at user-defined mode. \n"
		end
		if smooth != 0
			info *= "└─Lipschitz low-pass filter is applied. "
		end
		if ktol == 0
			info *= "No Krasny filter. "
		else
			info *= "Krasny filter with tolerance $ktol."
		end
		info *= "\nDiscretized with $(mesh.N) collocation points on [$(mesh.xmin), $(mesh.xmax)]."

		# Pre-allocate useful data
		k  = mesh.k
		x  = mesh.x
		∂ₓ	=  1im * k

		if dealias == 0
			Π⅔ 	= ones(size(k)) # no dealiasing (Π⅔=Id)
			Π 	= ones(size(k)) # no dealiasing (Π⅔=Id)
		else
			K = (mesh.kmax-mesh.kmin)/(3*dealias)
			Π⅔ 	= abs.(k) .<= K # Dealiasing low-pass filter
			if smooth ==0
				Π 	= abs.(k) .<= K # Dealiasing low-pass filter
			else
				Π  = max.( 0 , min.( 1 , 2/smooth*( 1 .-abs.(k)/K) )).^2	
			end
		end
		η = zeros(Float64, mesh.N)
        v = zeros(Float64, mesh.N)
		fftη = zeros(Complex{Float64}, mesh.N)
        fftv = zeros(Complex{Float64}, mesh.N)
		

		# Evolution equations are ∂t U = f(U)
		function f!(U)
			fftη .= U[1]
			fftv .= U[2]
			η .= real(ifft(fftη))
			v .= real(ifft(fftv))

			U[1] .= -∂ₓ.*( fftv .+ ϵ * Π.* fft( η .* v) )
			U[2] .= -∂ₓ.*( fftη .+ ϵ * Π.* fft( v.^2 )/2 )
			for u in U u[ abs.(u).< ktol ].=0 end
		end

		# Build raw data from physical data.
		# Discrete Fourier transform with, possibly, dealiasing and Krasny filter.
		function mapto(data::InitialData)
			U = [Π⅔.* fft(data.η(x)), Π⅔.*fft(data.v(x))]
			for u in U u[ abs.(u).< ktol ].=0 end
			return U
		end

		# Reconstruct physical variables from raw data
		# Return `(η,v,x)`, where
		# - `η` is the surface deformation;
		# - `v` is the derivative of the trace of the velocity potential;
		# - `x` is the vector of collocation points
		function mapfro(U)
			real(ifft(U[1])),real(ifft(U[2])),mesh.x
		end

		new(label, f!, mapto, mapfro, info)
    end
end

"""
	SaintVenant_fast(param;kwargs)

Same as `SaintVenant`, but faster.
"""
mutable struct SaintVenant_fast <: AbstractModel

	label 	:: String
	f!		:: Function
	mapto	:: Function
	mapfro	:: Function
	info	:: String

    function SaintVenant_fast(param::NamedTuple;
						mesh = Mesh(param),
						dealias=0,smooth=false,
						ktol=0,
						label="Saint-Venant"
						)

		# Set up
		ϵ 	= param.ϵ

		info = "Saint-Venant model.\n"
		info *= "├─Nonlinearity parameter ϵ=$(param.ϵ).\n"
		if dealias == 0
			info *= "└─No dealiasing. "
		elseif dealias == 1
			info *= "├─Dealiasing with Orszag's rule adapted to power 2 nonlinearity. \n"
		else
			info *= "├─Dealiasing at user-defined mode. \n"
		end
		if smooth != 0
			info *= "└─Lipschitz low-pass filter is applied. "
		end
		if ktol == 0
			info *= "No Krasny filter. "
		else
			@error("The Krasny filter will *not* be applied. Use 'SaintVenant' instead")
		end
		info *= "\nDiscretized with $(mesh.N) collocation points on [$(mesh.xmin), $(mesh.xmax)]."

		# Pre-allocate useful data
		k  = mesh.k
		x  = mesh.x
		∂ₓ	=  1im * k

		if dealias == 0
			Π⅔ 	= ones(size(k)) # no dealiasing (Π⅔=Id)
			Π 	= ones(size(k)) # no dealiasing (Π⅔=Id)
		else
			K = (mesh.kmax-mesh.kmin)/(3*dealias)
			Π⅔ 	= abs.(k) .<= K # Dealiasing low-pass filter
			if smooth ==0
				Π 	= abs.(k) .<= K # Dealiasing low-pass filter
			else
				Π  = max.( 0 , min.( 1 , 2/smooth*( 1 .-abs.(k)/K) )).^2	
			end
		end
		η = zeros(Complex{Float64}, mesh.N)
        v = zeros(Complex{Float64}, mesh.N)
		fftη = zeros(Complex{Float64}, mesh.N)
        fftv = zeros(Complex{Float64}, mesh.N)
		I = zeros(Complex{Float64}, mesh.N)

		Px  = plan_fft(η; flags = FFTW.MEASURE)
		ϵΠ=ϵ*Π
		∂=-∂ₓ

		# Evolution equations are ∂t U = f(U)
		function f!(U)
			fftη .= U[1]
			fftv .= U[2]

			ldiv!(η, Px, fftη )
			ldiv!(v, Px, fftv )


			η.*=v
			mul!(I, Px, η)
			I.*=ϵΠ
			fftv.+=I
			fftv.*=∂

			U[1] .= fftv

			v.*=v
			v./=2
			mul!(I, Px, v)
			I.*=ϵΠ
			fftη.+=I
			fftη.*=∂

			U[2] .= fftη
		end

		# Build raw data from physical data.
		# Discrete Fourier transform with, possibly, dealiasing and Krasny filter.
		function mapto(data::InitialData)
			U = [Π⅔.* fft(data.η(x)), Π⅔.*fft(data.v(x))]
			for u in U u[ abs.(u).< ktol ].=0 end
			return U
		end

		# Reconstruct physical variables from raw data
		# Return `(η,v,x)`, where
		# - `η` is the surface deformation;
		# - `v` is the derivative of the trace of the velocity potential;
		# - `x` is the vector of collocation points
		function mapfro(U)
			real(ifft(U[1])),real(ifft(U[2])),mesh.x
		end

		new(label, f!, mapto, mapfro, info)
    end
end
