export SaintVenant2D,SaintVenant2D_fast

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
- `dealias`: no dealisasing if set to `0` or `false` (default), otherwise `1/(3*dealias)` modes are set to `0` (corresponding to standard 2/3 Orszag rule if `dealias` is set to `1` or `true`);
- `smooth`: A smooth low-pass filter (whose scaling is defined by ) if set to `0` or `false` (default), otherwise only `2/(3*dealias)*(1-smooth/2)` modes are kept untouched;
- `label`: a label for future references (default is `"Saint-Venant"`);

# Return values
Generate necessary ingredients for solving an initial-value problem via `solve!`:
1. a function `SaintVenant2D.f!` to be called in explicit time-integration solvers;
2. a function `SaintVenant2D.mapto` which from `(η,vx,vy)` of type `InitialData` provides the raw data matrix on which computations are to be executed;
3. a function `SaintVenant2D.mapfro` which from such data matrix returns the Tuple of real vectors `(η,vx,vy,x,y)`, where
    - `η` is the values of surface deformation at collocation points `(x,y)`;
    - `vx,vy` are the velocity fields at collocation points `(x,y)`.

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
						smooth=false,
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
		kx=mesh.k;ky=mesh.k;x=mesh.x;y=mesh.x;nx=mesh.N;ny=mesh.N;
		∂y = 1im * ky' .* one.( x)
		∂x = one.(y') .* 1im .* kx
		
		if dealias == 0
			Π⅔ 	= ones(nx,ny) # no dealiasing (Π⅔=Id)
			Π 	= ones(nx,ny) # no dealiasing (Π⅔=Id)
		else
			K = (mesh.kmax-mesh.kmin)/(3*dealias)
			Π⅔ 	= (abs.(ky) .<= K)'.*(abs.(kx) .<= K) # Dealiasing low-pass filter
			if smooth ==0
				Px 	= abs.(kx) .<= K # Dealiasing low-pass filter
				Py 	= abs.(ky) .<= K # Dealiasing low-pass filter
			else
				Px  = max.( 0 , min.( 1 , 2/smooth*( 1 .-abs.(kx)/K) )).^2	
				Py  = max.( 0 , min.( 1 , 2/smooth*( 1 .-abs.(ky)/K) )).^2	
			end
			#P   = min.(( abs.(mesh.k) .<= K).*( (-abs.(mesh.k) .+ K)*dealias/mesh.dk ),1)	
			#P  =min.(( abs.(mesh.k) .<= K).*( (-abs.(mesh.k) .+ K)*dealias/mesh.dk .+1/2 ),1)	
			Π = Py'.*Px		
		end

		η = zeros(Float64, (nx,ny))
        vx = zeros(Float64, (nx,ny))
		vy = zeros(Float64, (nx,ny))
		fftη = zeros(Complex{Float64}, (nx,ny))
        fftvx = zeros(Complex{Float64}, (nx,ny))
		fftvy = zeros(Complex{Float64}, (nx,ny))
		
		# Evolution equations are ∂t U = f(U)
		function f!(U)
			fftη .= U[1]
			fftvx .= U[2]
			fftvy .= U[3]

			η .= real(ifft(fftη))
			vx .= real(ifft(fftvx))
			vy .= real(ifft(fftvy))

			U[1] .= -∂x.*( fftvx .+ ϵ * Π.* fft( η .* vx) ) +
						-∂y.*( fftvy .+ ϵ * Π.* fft( η .* vy) )
			if hamiltonian
				U[2] .= -∂x.*( fftη .+ ϵ/2 * Π.* fft( vx.^2+vy.^2 ) )
				U[3] .= -∂y.*( fftη .+ ϵ/2 * Π.* fft( vx.^2+vy.^2 ) )
			else
				U[2] .= -∂x.*( fftη ) .- ϵ * Π.* fft( vx.* ifft(∂x.*fftvx) + vy.* ifft(∂y.*fftvx) )
				U[3] .= -∂y.*( fftη ) .- ϵ * Π.* fft( vx.* ifft(∂x.*fftvy) + vy.* ifft(∂y.*fftvy) )
			end		
			for u in U u[ abs.(u).< ktol ].=0 end
		end

		# Build raw data from physical data.
		# Discrete Fourier transform with, possibly, dealiasing and Krasny filter.
		function mapto(data::InitialData)
			#U = cat(Π⅔.* fft(data.η(x,y)), Π⅔.*fft(data.vx(x,y)), Π⅔.*fft(data.vy(x,y)),dims=3)
			U = [Π⅔.* fft(data.η(x,y)), Π⅔.*fft(data.vx(x,y)), Π⅔.*fft(data.vy(x,y))]

			for u in U u[ abs.(u).< ktol ].=0 end
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

"""
	SaintVenant_fast(param;kwargs)

Same as `SaintVenant`, but faster.
"""
mutable struct SaintVenant2D_fast <: AbstractModel

	label 	:: String
	f!		:: Function
	mapto	:: Function
	mapfro	:: Function
	info	:: String

    function SaintVenant2D_fast(param::NamedTuple;
						mesh = Mesh(param),
						dealias=0,ktol=0,
						smooth=false,
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
		kx=mesh.k;ky=mesh.k;x=mesh.x;y=mesh.x;nx=mesh.N;ny=mesh.N;
		
		if dealias == 0
			Π⅔ 	= ones(nx,ny) # no dealiasing (Π⅔=Id)
			Π 	= ones(nx,ny) # no dealiasing (Π⅔=Id)
		else
			K = (mesh.kmax-mesh.kmin)/(3*dealias)
			Π⅔ 	= (abs.(ky) .<= K)'.*(abs.(kx) .<= K) # Dealiasing low-pass filter
			if smooth ==0
				Px 	= abs.(kx) .<= K # Dealiasing low-pass filter
				Py 	= abs.(ky) .<= K # Dealiasing low-pass filter
			else
				Px  = max.( 0 , min.( 1 , 2/smooth*( 1 .-abs.(kx)/K) )).^2	
				Py  = max.( 0 , min.( 1 , 2/smooth*( 1 .-abs.(ky)/K) )).^2	
			end
			#P   = min.(( abs.(mesh.k) .<= K).*( (-abs.(mesh.k) .+ K)*dealias/mesh.dk ),1)	
			#P  =min.(( abs.(mesh.k) .<= K).*( (-abs.(mesh.k) .+ K)*dealias/mesh.dk .+1/2 ),1)	
			Π = Py'.*Px		
		end
		
    	f  = zeros(ComplexF64,(nx,ny))
    	f̂  = similar(f)
    	fᵗ = zeros(ComplexF64,(ny,nx))
    	f̂ᵗ = similar(fᵗ)

    
    	FFTW.set_num_threads(4)
    	Px = plan_fft(f,  1, flags=FFTW.PATIENT)    
    	Py = plan_fft(fᵗ, 1, flags=FFTW.PATIENT)
		P2 = plan_fft(f,  2, flags=FFTW.PATIENT)    

    
		η = zeros(Float64, (nx,ny))
        vx = zeros(Float64, (nx,ny))
		vy = zeros(Float64, (nx,ny))
		fftη = zeros(Complex{Float64}, (nx,ny))
        fftvx = zeros(Complex{Float64}, (nx,ny))
		fftvy = zeros(Complex{Float64}, (nx,ny))
		I = zeros(Complex{Float64}, mesh.N)

		ϵΠ=ϵ*Π
		∂=-∂ₓ

		function my_ifft!(fft,f,fᵗ,gᵗ)
			ldiv!(f, Px, fft )
			transpose!(fᵗ,f)
			ldiv!(gᵗ, Py, fᵗ )
			transpose!(fft,gᵗ)
		end

		function my_ifft2!(fft,f)
			ldiv!(f, Px, fft )
			ldiv!(fft, P2, f )
		end



		# Evolution equations are ∂t U = f(U)
		function f!(U)
			fftη .= U[1]
			fftvx .= U[2]
			fftvy .= U[3]

			# TO DO Complete Saint-Venant using
			# https://pnavaro.github.io/math-julia/diffeq/rotation-with-fft.html

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
