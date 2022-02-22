export Matsuno_fast,Matsuno

"""
	Matsuno_fast(param;dealias,label)

Same as `Matsuno`, but faster.
"""
mutable struct Matsuno_fast <: AbstractModel

	label   :: String
	f!		:: Function
	mapto	:: Function
	mapfro	:: Function
	info	:: String

    function Matsuno_fast(param::NamedTuple;
							IL	    = false,
							dealias=false,
							label = "Matsuno" )

		# Set-up
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
			IL = true;  # IL (=Infinite layer) is a flag to be used thereafter
			μ = 1; ν = 1; # Then we should set μ=ν=1 in subsequent formula.
		end
		mesh = Mesh(param)

		# Print information
		info = "Matsuno model.\n"
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
		info *= "\nDiscretized with $(mesh.N) collocation points on [$(mesh.xmin), $(mesh.xmax)]."

		# Pre-allocate useful data
		x = mesh.x
		k = mesh.k
		if IL == true
			Tμ 	= -1im * sign.(k)
			G₀ 	= abs.(k)
		else
			Tμ 	= -1im* sign.(k) .* tanh.(sqrt(μ)*abs.(k))
			G₀ 	= sqrt(μ)*abs.(k).*tanh.(sqrt(μ)*abs.(k))
		end
		∂ₓ	=  1im * sqrt(μ)* k            # Differentiation
		if dealias == true || dealias == 1
			Π⅔    = abs.(k) .< (mesh.kmax-mesh.kmin)/3    # Dealiasing low-pass filter
		else
			Π⅔    = zero(k) .+ 1     		# No dealiasing (Π⅔=Id)
		end

        hnew = zeros(Complex{Float64}, mesh.N)
        unew = zeros(Complex{Float64}, mesh.N)

        I₀ = zeros(Complex{Float64}, mesh.N)
        I₁ = zeros(Complex{Float64}, mesh.N)
        I₂ = zeros(Complex{Float64}, mesh.N)
        I₃ = zeros(Complex{Float64}, mesh.N)

        Px  = plan_fft(hnew; flags = FFTW.MEASURE)

		# Evolution equations are ∂t U = f(U)
		function f!(U)

		    for i in eachindex(hnew)
		        hnew[i] = G₀[i] * U[i,1]
		    end

		    ldiv!(unew, Px, hnew )

		    for i in eachindex(hnew)
		        hnew[i] = ∂ₓ[i] * U[i,1]
		    end

		    ldiv!(I₁, Px, hnew)

		    unew  .*= I₁

		    mul!(I₁, Px, unew)

		    I₁  .*= ϵ .* Π⅔
		    I₁  .-= hnew

		    ldiv!(hnew, Px, view(U,:,1))
		    ldiv!(unew, Px, view(U,:,2))

		    I₂    .= hnew .* unew

		    mul!(I₃, Px, I₂)

		    I₃    .*= ∂ₓ

		    for i in eachindex(Tμ)
		        U[i,1]  = Tμ[i] * U[i,2]
		        I₀[i] = G₀[i] * U[i,2]
		    end

		    ldiv!(I₂, Px, I₀)

		    I₂    .*= hnew

		    mul!(hnew, Px, I₂)

		    hnew  .*= Tμ
		    I₃    .+= hnew
		    I₃    .*= ϵ .* Π⅔

		    for i in eachindex(I₃)
		        U[i,1] -= I₃[i]
		    end
			U[:,1] ./= sqrt(μ)/ν

		    I₃    .=  unew.^2

		    mul!(unew, Px, I₃)

		    unew  .*= ∂ₓ
		    unew  .*= ϵ/2/ν .* Π⅔
		    I₁    .-= unew

		    for i in eachindex(I₁)
		        U[i,2] =  I₁[i]/sqrt(μ)
		    end

		end

		# Build raw data from physical data.
		function mapto(data::InitialData)
			fftη = Π⅔ .* fft(data.η(x));
			fftv = Π⅔ .* fft(data.v(x));
			U = [fftη fftv-ϵ* Π⅔ .*fft(ifft(Tμ.*fftv).*ifft(∂ₓ.*fftη) )]
			return U
		end

		# Take raw data and return `(η,v)`, where
		# - `η` is the surface deformation;
		# - `v` is the velocity variable.
		function mapfro(U)
			#real(ifft(view(U,:,1))),real(ifft(view(U,:,2))+ϵ* ifft(Tμ.*view(U,:,2)).*ifft(∂ₓ.*view(U,:,1)))
			real(ifft(U[:,1])),real(ifft(U[:,2])+ϵ* ifft(Tμ.*U[:,2]).*ifft(∂ₓ.*U[:,1]))

		end

		new(label, f!, mapto, mapfro, info )
    end
end

"""
	Matsuno(param;kwargs)

Define an object of type `AbstractModel` in view of solving the initial-value problem for
the quadratic deep-water model proposed by [Matsuno](https://doi.org/10.1103/PhysRevLett.69.609).

# Argument
`param` is of type `NamedTuple` and must contain
- dimensionless parameters `ϵ` (nonlinearity) and `μ` (dispersion);
- optionally, `ν` the shallow/deep water scaling factor. By default, `ν=1` if `μ≦1` and `ν=1/√μ` otherwise. Set the infinite-layer case if `ν=0`, or `μ=Inf`.
- numerical parameters to construct the mesh of collocation points as `mesh = Mesh(param)`

## Optional keyword arguments
- `IL`: Set the infinite-layer case if `IL=true` (or `μ=Inf`, or `ν=0`), in which case `ϵ` is the steepness parameter. Default is `false`.
- `dealias`: dealiasing with `1/3` Orlicz rule if `true` or no dealiasing if `false` (by default);
- `label`: a label for future references (default is `"Matsuno"`);

# Return values
Generate necessary ingredients for solving an initial-value problem via `solve!`:
1. a function `DeepQuadratic.f!` to be called in explicit time-integration solvers;
2. a function `DeepQuadratic.mapto` which from `(η,v)` of type `InitialData` provides the raw data matrix on which computations are to be executed;
3. a function `DeepQuadratic.mapfro` which from such data matrix returns the Tuple of real vectors `(η,v)`, where
    - `η` is the surface deformation;
    - `v` is a velocity variable which is *not* the derivative of the trace of the velocity potential (if not null).

"""
mutable struct Matsuno <: AbstractModel

	label   :: String
	f!		:: Function
	mapto	:: Function
	mapfro	:: Function
	info 	:: String


    function Matsuno(param::NamedTuple ;
						IL	    = false,
						dealias=false,
						label = "Matsuno" )

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
			IL = true;  # IL (=Infinite layer) is a flag to be used thereafter
			μ = 1; ν = 1; # Then we should set μ=ν=1 in subsequent formula.
		end
		mesh = Mesh(param)

		# Print information
		info = "Matsuno model.\n"
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
		info *= "\nDiscretized with $(mesh.N) collocation points on [$(mesh.xmin), $(mesh.xmax)]."

		# Pre-allocate useful data
		x = mesh.x
		k = mesh.k
		if IL == true
			Tμ 	= -1im * sign.(k)
			G₀ 	= abs.(k)
		else
			Tμ 	= -1im* sign.(k) .* tanh.(sqrt(μ)*abs.(k))
			G₀ 	= sqrt(μ)*abs.(k).*tanh.(sqrt(μ)*abs.(k))
		end
		∂ₓ	=  1im * sqrt(μ)* k            # Differentiation

		if dealias == true || dealias == 1
			Π⅔    = abs.(k) .< (mesh.kmax-mesh.kmin)/3 	# Dealiasing low-pass filter
		else
			Π⅔    = zero(k) .+ 1     		# No dealiasing (Π⅔=Id)
		end

        hnew = zeros(Complex{Float64}, mesh.N)
        unew = zeros(Complex{Float64}, mesh.N)

        I₁ = zeros(Complex{Float64}, mesh.N)
        I₂ = zeros(Complex{Float64}, mesh.N)
        I₃ = zeros(Complex{Float64}, mesh.N)

		# Evolution equations are ∂t U = f(U)
		function f!(U)

		   hnew .= ifft(U[:,1])
		   unew .= ifft(U[:,2])
		   I₃ .= fft(ifft(∂ₓ.*U[:,1]).*ifft(G₀.*U[:,1]))
		   I₁ .= Tμ.*U[:,2].-ϵ*Π⅔.*(Tμ.*fft(hnew.*ifft(G₀.*U[:,2])).+∂ₓ.*fft(hnew.*unew))
		   I₂ .= -ν*(∂ₓ.*U[:,1])+ϵ*Π⅔.*(ν*I₃-∂ₓ.*fft(unew.^2)/2)
		   #
		   U[:,1] .= I₁/sqrt(μ)/ν
		   U[:,2] .= I₂/sqrt(μ)/ν

		end

		# Build raw data from physical data.
		function mapto(data::InitialData)
			fftη = Π⅔ .* fft(data.η(x));
			fftv = Π⅔ .* fft(data.v(x));
			U = [fftη fftv-ϵ* Π⅔ .*fft(ifft(Tμ.*fftv).*ifft(∂ₓ.*fftη) )]
			return U
		end

		# Take raw data and return `(η,v)`, where
		# - `η` is the surface deformation;
		# - `v` is the velocity variable.
		function mapfro(U)
			real(ifft(U[:,1])),real(ifft(U[:,2])+ϵ* ifft(Tμ.*U[:,2]).*ifft(∂ₓ.*U[:,1]))
		end

		new(label, f!, mapto, mapfro, info )
    end
end
