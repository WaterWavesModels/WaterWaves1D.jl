export modifiedMatsuno

"""
    modifiedMatsuno(param;kwargs)

Define an object of type `AbstractModel` in view of solving the initial-value problem for
the modified Matsuno model

# Argument
`param` is of type `NamedTuple` and must contain
- dimensionless parameters `ϵ` (nonlinearity) and `μ` (dispersion);
- optionally, `ν` the shallow/deep water scaling factor. By default, `ν=1` if `μ≦1` and `ν=1/√μ` otherwise. Set the infinite-layer case if `ν=0`, or `μ=Inf`.
- numerical parameters to construct the mesh of collocation points, if `mesh` is not provided as a keyword argument.

## Optional keyword arguments
- `IL`: Set the infinite-layer case if `IL=true` (or `μ=Inf`, or `ν=0`), in which case `ϵ` is the steepness parameter. Default is `false`.
- `mesh`: the mesh of collocation points. By default, `mesh = Mesh(param)`;
- `ktol`: tolerance of the low-pass Krasny filter (default is `0`, i.e. no filtering);
- `dealias`: dealiasing with Orlicz rule `1-dealias/(dealias+2)` (default is `0`, i.e. no dealiasing);
- `label`: a label for future references (default is `"modified Matsuno"`);

# Return values
Generate necessary ingredients for solving an initial-value problem via `solve!`:
1. a function `modifiedMatsuno.f!` to be called in explicit time-integration solvers;
2. a function `modifiedMatsuno.mapto` which from `(η,v)` of type `InitialData` provides the raw data matrix on which computations are to be executed;
3. a function `modifiedMatsuno.mapfro` which from such data matrix returns the Tuple of real vectors `(η,v,x)`, where
    - `η` is the values of surface deformation at collocation points `x`;
    - `v` is the derivative of the trace of the velocity potential at `x`.

"""
mutable struct modifiedMatsuno <: AbstractModel

	label   :: String
	f!		:: Function
	mapto	:: Function
	mapfro	:: Function
	info	:: String

    function modifiedMatsuno(param::NamedTuple;
							mesh = Mesh(param),
							IL	    = false,
							ktol	= 0,
							dealias	= 0,
							label	= "modified Matsuno"
							)

		# Set up
		μ 	= param.μ
		ϵ 	= param.ϵ
		if !in(:ν,keys(param))
			if μ > 1
				ν = 1/sqrt(μ)
				nu = "1/√μ (deep water case)"
			else
				ν = 1
				nu = "1 (shallow water case)"
			end
		else
			ν = param.ν
			nu = "$ν"
		end
		if μ == Inf || ν==0 || IL == true # infinite layer case
			IL = true;  # IL (=Infinite layer) is a flag to be used therefafter
			μ = 1; ν = 1; # Then we should set μ=ν=1 in subsequent formula.
		else IL = false
		end

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
			Tμ 	= -1im * sign.(k)
			G₀ 	= abs.(k)
		else
        	Tμ 	= -1im* sign.(k) .* tanh.(sqrt(μ)*abs.(k))
			G₀ 	= sqrt(μ)*abs.(k).*tanh.(sqrt(μ)*abs.(k))
		end
		∂ₓ	=  1im * sqrt(μ)* k            # Differentiation
		z = zeros(Complex{Float64}, mesh.N)
		η = copy(z) ; v = copy(z) ; F = copy(z) ;
		fftη = copy(z) ; fftv = copy(z) ; Q = copy(z) ; R = copy(z) ;

		# Build raw data from physical data.
		function mapto(data::InitialData)
			fftη .= Π⅔ .* fft(data.η(x));
			fftv .= Π⅔ .* fft(data.v(x));
			U = [fftη fftv-ϵ* Π⅔ .*fft(ifft(Tμ.*fftv).*ifft(∂ₓ.*fftη) )]
			U[abs.(U).< ktol ].=0
			return U
		end

		# Reconstruct physical variables from raw data
		# Return `(η,v,x)`, where
		# - `η` is the surface deformation;
		# - `v` is the derivative of the trace of the velocity potential;
		# - `x` is the vector of collocation points
		function mapfro(U;n=10)
			∂ζ=ifft(∂ₓ.*U[:,1]);
			z.=U[:,2];v.=U[:,2];
			for j=1:n
				v.=z+ϵ*Π⅔ .* fft( ∂ζ .* ifft(Tμ.*v))
			end
			real(ifft(U[:,1])),real(ifft(v)),mesh.x
		end

		# Evolution equations are ∂t U = f(U)
		function f!(U)

			fftη .= U[:,1]
		    fftv .= U[:,2]
			η  .= ifft(U[:,1])
			v  .= ifft(U[:,2])
			Q .= Tμ.*fftv
			Q -= ϵ *∂ₓ.*fft(η.*ifft(fftv)) .+ ϵ*Tμ.*fft(η.*ifft(G₀.*fftv))
			F .= exp.(-ϵ*ifft(G₀.*fftη))
			R .= -fft(F.*ifft(∂ₓ.*fftη))*ν
			R -= ϵ/2*∂ₓ.*fft(v.^2)

			U[:,1] .= Π⅔.*Q/sqrt(μ)/ν
		   	U[:,2] .= Π⅔.*R/sqrt(μ)/ν
			U[abs.(U).< ktol ].=0

		end


        new(label, f!, mapto, mapfro, info )
    end
end
