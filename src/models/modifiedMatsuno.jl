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
- `ktol`: tolerance of the low-pass Krasny filter (default is `0`, i.e. no filtering);
- `dealias`: dealiasing with Orlicz rule `1-dealias/(dealias+2)` (default is `0`, i.e. no dealiasing);
- `label`: a label for future references (default is `"modified Matsuno"`);
- `verbose`: prints information if `true` (default is `true`).

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

    function modifiedMatsuno(param::NamedTuple;
							ν		=nothing,
							ktol	=0,
							dealias	=0,
							label	="modified Matsuno",
							verbose	=true )

		if verbose @info "Build the modified Matsuno model." end

		μ 	= param.μ
		ϵ 	= param.ϵ
		if ν == nothing
			if μ > 1
				ν = 1/sqrt(μ)
			else
				ν = 1
			end
		end

		mesh = Mesh(param)
		x = mesh.x
		k = copy(mesh.k)
		K = mesh.kmax * (1-dealias/(2+dealias))
		Π⅔ 	= abs.(mesh.k) .<= K # Dealiasing low-pass filter
		if dealias == 0
			if verbose @info "no dealiasing" end
			Π⅔ 	= ones(size(mesh.k))
		elseif verbose
			@info "dealiasing : spectral scheme for power  $(dealias + 1) nonlinearity"
		end
		if μ == Inf || ν==0
			∂ₓF₀ 	= 1im * sign.(mesh.k)
			G₀ 	= abs.(mesh.k)
			μ = 1
			ν = 1
		else
        	∂ₓF₀ 	= 1im* sign.(mesh.k) .* tanh.(sqrt(μ)*abs.(mesh.k))
			G₀ 	= sqrt(μ)*abs.(mesh.k).*tanh.(sqrt(μ)*abs.(mesh.k))
		end
		∂ₓ	=  1im * sqrt(μ)* mesh.k            # Differentiation
		z = zeros(Complex{Float64}, mesh.N)
		η = copy(z) ; v = copy(z) ; F = copy(z) ;
		fftη = copy(z) ; fftv = copy(z) ; Q = copy(z) ; R = copy(z) ;

		# Discrete Fourier transform with, possibly, dealiasing and Krasny filter.
		function mapto(data::InitialData)
			U = [Π⅔ .* fft(data.η(x)) Π⅔ .*fft(data.v(x))]
			U[abs.(U).< ktol ].=0
			return U
		end

		# Return `(η,v)`, where
		# - `η` is the surface deformation;
		# - `v` is the derivative of the trace of the velocity potential.
		# Inverse Fourier transform and takes the real part.
		function mapfro(U)
			real(ifft(U[:,1])),real(ifft(U[:,2]))
		end

		# Evolution equations are ∂t U = f!(U)
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


        new(label, f!, mapto, mapfro )
    end
end
