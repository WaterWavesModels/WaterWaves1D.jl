export SVN

"""
    SVN(param;kwargs)

Define an object of type `AbstractModel` in view of solving the initial-value problem for
a the multilayer Saint-Venant model, possibly with viscosity and/or diffusivity.

# Argument
`param` is of type `NamedTuple` and must contain
- 'M': the number of layers in the model;
- dimensionless parameters `ϵ` (nonlinearity), `ν` (viscosity) and `κ` (diffusivity);
- 'ρ' a `Vector` of size `M`, the density of each layer from top to bottom;
- 'δ' a `Vector` of size `M`, the depth of each layer from top to bottom;
- 'u₀' a `Vector` of size `M`, the horizontal velocity of the background shear flow at each layer from top to bottom;
- numerical parameters to construct the mesh of collocation points, if `mesh` is not provided as a keyword argument.

## Optional keyword arguments
- `mesh`: the mesh of collocation points. By default, `mesh = Mesh(param)`;
- `ktol`: tolerance of the low-pass Krasny filter (default is `0`, i.e. no filtering);
- `dealias`: dealiasing with Orlicz rule `1-dealias/(dealias+2)` (default is `0`, i.e. no dealiasing);
- `label`: a label for future references (default is `"SVN"`);


# Return values
Generate necessary ingredients for solving an initial-value problem via `solve!`:
1. a function `WhithamBoussinesq.f!` to be called in explicit time-integration solvers;
2. a function `WhithamBoussinesq.mapto` which from `(η,v)` of type `InitialData` provides the raw data matrix on which computations are to be executed.
3. a function `WhithamBoussinesq.mapfro` which from such data matrix returns the Tuple of real vectors `(η,v,x)`, where
	- `η` is the values of surface deformation at collocation points `x`;
	- `v` is the derivative of the trace of the velocity potential at `x`.

"""
mutable struct SVN <: AbstractModel

	label   :: String
	f!		:: Function
	mapto	:: Function
	mapfro	:: Function
	info 	:: String

    function SVN(param::NamedTuple;
								mesh = Mesh(param),
								dealias = 0,
								ktol	= 0,
								label 	= nothing
								)

		# Set up
		κ 	= param.κ
		ϵ 	= param.ϵ
		ρ	= param.ρ
		δ	= param.δ
		ν	= param.ν

		u₀	 = param.u₀

		M=length(ρ)
		N=mesh.N

		if isnothing(label) label = "SVN" end
		
		# Print information
		info = "$label model with N=$M layers.\n"
		info *= "├─Diffusivity parameter κ=$κ, viscosity parameter ν=$ν, nonlinearity parameter ϵ=$ϵ.\n"
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
		
		if dealias == 0
			Π⅔ 	= ones(size(k)) # no dealiasing (Π⅔=Id)
		else
			K = (mesh.kmax-mesh.kmin)/(2+dealias)
			Π⅔ 	= abs.(k) .<= K # Dealiasing low-pass filter
		end
		η = zeros(Float64, mesh.N,M)
        v = zeros(Float64, mesh.N,M)
		fftη = zeros(Complex{Float64}, mesh.N,M)
        fftv = zeros(Complex{Float64}, mesh.N,M)


		# Evolution equations are ∂t U = f(U)
		function f!(U)

			fftη .= U[:,1:M]
			fftv .= U[:,M+1:2*M]

			for m = 1:M
				v[:,m] .= real(ifft(fftv[:,m]))
				η[:,m] .= real(ifft(fftη[:,m]))
			end
			for m = 1:M

				U[:,m  ] .= -∂ₓ.*(u₀[m].*fftη[:,m].+δ[m].*fftv[:,m].+ϵ*Π⅔.*fft(η[:,m].*v[:,m]))+κ*∂ₓ.*∂ₓ.*fftη[:,m]
		   		U[:,m+M] .= -∂ₓ.*(u₀[m].*fftv[:,m].+ϵ/2*Π⅔.*fft(v[:,m].^2))+ν*∂ₓ.*∂ₓ.*fftv[:,m]
				for n = 1:M 
					U[:,m+M] .-= 1/M*min(ρ[n]/ρ[m],1)*∂ₓ.*fftη[:,n]
				end
			end
			U[abs.(U).< ktol ].=0

		end

		# Build raw data from physical data.
		# Discrete Fourier transform with, possibly, dealiasing and Krasny filter.
		function mapto(data::InitialData)
			U = zeros(Complex{Float64}, mesh.N,2*M)

			for m = 1:M
				U[:,m  ] .= Π⅔ .* fft(data.η(x)[:,m]) 
				U[:,m+M] .= Π⅔ .* fft(data.v(x)[:,m])
			end
			U[abs.(U).< ktol ].=0
			return U
		end

		# Reconstruct physical variables from raw data
		# Return `(η,v,x)`, where
		# - `η` is the vector of layer depth deviations;
		# - `v` is the vector of horizontal velocities;
		# - `x` is the vector of collocation points
		function mapfro(U)
			η = zeros(Float64, mesh.N,M)
        	v = zeros(Float64, mesh.N,M)
			for m = 1:M
				η[:,m]=real(ifft(U[:,m]))
				v[:,m]=real(ifft(U[:,m+M]))
			end
			return η,v,mesh.x
		end

        new(label, f!, mapto, mapfro, info )
    end
end
