export modifiedMatsuno

"""
    modifiedMatsuno(param;kwargs)

Define an object of type `AbstractModel` in view of solving the initial-value problem for
the modified Matsuno model

# Argument
`param` is of type `NamedTuple` and must contain
- dimensionless parameters `ϵ` (nonlinearity) and `μ` (dispersion);
- numerical parameters to construct the mesh of collocation points as `mesh = Mesh(param)`

## Optional keyword arguments
- `ν`: shallow/deep water multiplication factor. By default, `ν=1` if `μ≦1` and `ν=1/√μ` otherwise. Set the infinite-layer case if `ν=0` (or `μ=Inf`).
- `IL`: Set the infinite-layer case if `IL=true` (or `μ=Inf`, or `ν=0`), in which case `ϵ` is the steepness parameter. Default is `false`.
- `ktol`: tolerance of the low-pass Krasny filter (default is `0`, i.e. no filtering);
- `dealias`: dealiasing with Orlicz rule `1-dealias/(dealias+2)` (default is `0`, i.e. no dealiasing);
- `label`: a label for future references (default is `"modified Matsuno"`);

# Return values
Generate necessary ingredients for solving an initial-value problem via `solve!`:
1. a function `modifiedMatsuno.f!` to be called in explicit time-integration solvers;
2. a function `modifiedMatsuno.mapto` which from `(η,v)` of type `InitialData` provides the raw data matrix on which computations are to be executed;
3. a function `modifiedMatsuno.mapfro` which from such data matrix returns the Tuple of real vectors `(η,v)`, where
    - `η` is the surface deformation;
    - `v` is the derivative of the trace of the velocity potential.

"""
mutable struct modifiedMatsuno <: AbstractModel

	label   :: String
	f!		:: Function
	mapto	:: Function
	mapfro	:: Function
	info	:: String

    function modifiedMatsuno(param::NamedTuple;
							ν		= nothing,
							IL	    = false,
							ktol	= 0,
							dealias	= 0,
							label	= "modified Matsuno"
							)

		# Set up
		μ 	= param.μ
		ϵ 	= param.ϵ
		if isnothing(ν)
			if μ > 1
				ν = 1/sqrt(μ)
				nu = "1/√μ (deep water case)"
			else
				ν = 1
				nu = "1 (shallow water case)"
			end
		else
			nu = "$ν"
		end
		if μ == Inf || ν==0 # infinite layer case
			IL = true;  # IL (=Infinite layer) is a flag to be used therefafter
			μ = 1; ν = 1; # Then we should set μ=ν=1 in subsequent formula.
		else IL = false
		end
		mesh = Mesh(param)

		# Print information
		info = "Modified 'exponential' Matsuno model of Duchêne and Melinand.\n"
		if IL == true
			info *= "├─Steepness parameter ϵ=$ϵ (infinite depth case).\n"
		else
			info *= "├─Shallowness parameter μ=$μ, nonlinearity parameter ϵ=$ϵ, \
					scaling parameter ν=$nu.\n"
		end
		if dealias == true || dealias == 1
			info *= "└─Dealiasing with Orszag’s 3/2 rule. "
		else
			info *= "└─No dealiasing. "
		end
		if ktol == 0
			info *= "No Krasny filter. "
		else
			info *= "Krasny filter with tolerance $ktol."
		end
		info *= "\nDiscretized with $(mesh.N) collocation points on [$(mesh.xmin), $(mesh.xmax)]."
		@warn "The velocity is consistent with the \
		derivative of the trace of the velocity potential \
		used for water waves only when they are null."

		# Pre-allocate useful data
		x = mesh.x
		k = mesh.k
		if dealias == 0
			Π⅔ 	= ones(size(k)) # no dealiasing (Π⅔=Id)
		else
			K = (mesh.kmax-mesh.kmin)/(2+dealias)
			Π⅔ 	= abs.(k) .<= K # Dealiasing low-pass filter
		end
		if IL == true
			∂ₓF₀ 	= 1im * sign.(k)
			G₀ 	= abs.(k)
		else
        	∂ₓF₀ 	= 1im* sign.(k) .* tanh.(sqrt(μ)*abs.(k))
			G₀ 	= sqrt(μ)*abs.(k).*tanh.(sqrt(μ)*abs.(k))
		end
		∂ₓ	=  1im * sqrt(μ)* k            # Differentiation
		z = zeros(Complex{Float64}, mesh.N)
		η = copy(z) ; v = copy(z) ; F = copy(z) ;
		fftη = copy(z) ; fftv = copy(z) ; Q = copy(z) ; R = copy(z) ;

		# Build raw data from physical data.
		# Discrete Fourier transform with, possibly, dealiasing and Krasny filter.
		# This function is correct only when the initial data for v is zero.
		function mapto(data::InitialData)
			U = [Π⅔ .* fft(data.η(x)) Π⅔ .*fft(data.v(x))]
			U[abs.(U).< ktol ].=0
			return U
		end

		# Return `(η,v)`, where
		# - `η` is the surface deformation;
		# - `v` is the velocity variable.
		# Inverse Fourier transform and takes the real part.
		function mapfro(U)
			real(ifft(U[:,1])),real(ifft(U[:,2]))
		end

		# Evolution equations are ∂t U = f(U)
		function f!(U)

			fftη .= U[:,1]
		    fftv .= U[:,2]
			η  .= ifft(U[:,1])
			v  .= ifft(U[:,2])
			Q .= -∂ₓF₀.*fftv
			Q += -ϵ *∂ₓ.*fft(η.*ifft(fftv)) .+ ϵ*∂ₓF₀.*fft(η.*ifft(G₀.*fftv))
			F .= exp.(-ϵ*ifft(G₀.*fftη))
			R .= -fft(F.*ifft(∂ₓ.*fftη))*ν
			R += ϵ/2*∂ₓ.*fft(-v.^2)

			U[:,1] .= Π⅔.*Q/sqrt(μ)/ν
		   	U[:,2] .= Π⅔.*R/sqrt(μ)/ν
			U[abs.(U).< ktol ].=0

		end


        new(label, f!, mapto, mapfro, info )
    end
end
