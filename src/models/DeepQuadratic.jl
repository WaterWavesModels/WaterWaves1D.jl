export DeepQuadratic_fast,DeepQuadratic

"""
    DeepQuadratic_fast(param;dealias,label,verbose)

Same as `DeepQuadratic`, but faster.
"""
mutable struct DeepQuadratic_fast <: AbstractModel

	label   :: String
	f!		:: Function
	mapto	:: Function
	mapfro	:: Function

    function DeepQuadratic_fast( param::NamedTuple;
						dealias=false, label="deep quadratic", verbose=true )

		# Set up
		ϵ = param.ϵ

		if verbose  # Print information
			info = "Build the deep quadratic model.\n"
			info *= "Steepness parameter ϵ=$ϵ (infinite depth case).\n"
			if dealias == true || dealias == 1
				info *= "Dealiasing with Orszag’s 3/2 rule. "
			else
				info *= "No dealiasing. "
			end
			info *= "\nShut me up with keyword argument `verbose = false`."
			@info info
			@warn "You should provide the initial value for v by ∂t η = - ∂x v.\n\
		It equals the derivative of the trace of the velocity potential \
		used for water waves only when they are null."
		end

		# Pre-allocate useful data
		mesh  = Mesh(param)
        x = mesh.x
		k = mesh.k
        Γ = abs.(k)
        Dx    =  1im * k            	# Differentiation
        H     = -1im * sign.(k)     	# Hilbert transform
		if dealias == true || dealias == 1
			Π⅔    = Γ .< (mesh.kmax-mesh.kmin)/3 	# Dealiasing low-pass filter
		else
			Π⅔    = zero(Γ) .+ 1     	 	# No dealiasing (Π⅔=Id)
		end

        hnew = zeros(Complex{Float64}, mesh.N)
        unew = zeros(Complex{Float64}, mesh.N)

        I₁ = zeros(Complex{Float64}, mesh.N)
        I₂ = zeros(Complex{Float64}, mesh.N)
        I₃ = zeros(Complex{Float64}, mesh.N)

        Px  = plan_fft(hnew; flags = FFTW.MEASURE)

		# Evolution equations are ∂t U = f(U)
		function f!(U)

		    ldiv!(hnew, Px , view(U,:,1))

		    I₁  .= view(U,:,2) .* Γ
		    ldiv!(unew, Px , I₁)
		    unew .^= 2
		    mul!(I₁, Px , unew)
		    I₁ .*= H

		    I₂  .= view(U,:,1) .* Dx
		    ldiv!(unew, Px , I₂)
		    unew .*= hnew
		    mul!(I₂, Px , unew)

		    I₃  .= view(U,:,1) .* Γ
		    ldiv!(unew, Px, I₃)
		    unew .*= hnew
		    mul!(I₃ , Px , unew)
		    I₃ .*= H

		    hnew  .= .- view(U,:,2) .* Dx

		    I₁ .-= I₂
		    I₁ .-= I₃
		    I₁ .*= Π⅔
		    I₁ .*= ϵ

		    U[:,2] .= view(U,:,1) .* H .+ I₁
		    U[:,1] .= hnew

		end

		# Build raw data from physical data.
		# Discrete Fourier transform with, possibly, dealiasing.
		# This function is correct only when the initial data for v is zero.
		function mapto(data::InitialData)
		    [Π⅔ .* fft(data.η(x)) Π⅔ .*fft(data.v(x))]
		end

		# Return `(η,v)`, where
		# - `η` is the surface deformation;
		# - `v` is detemined by ∂t η + ∂x v = 0.
		# Inverse Fourier transform and takes the real part.
		function mapfro(U)
			real(ifft(U[:,1])),real(ifft(U[:,2]))
		end
		new(label, f!, mapto, mapfro )

    end
end

"""
    DeepQuadratic(param;dealias,label,verbose)

Define an object of type `AbstractModel` in view of solving the initial-value problem for
the quadratic deep-water model proposed by [Akers and Milewski](https://doi.org/10.1137/090758386)
and [Cheng, Granero-Belinchón, Shkoller and Milewski](https://doi.org/10.1007/s42286-019-00005-w)

# Arguments
`param` is of type `NamedTuple` and must contain
- the dimensionless parameters `ϵ` (nonlinearity);
- numerical parameters to construct the mesh of collocation points as `mesh = Mesh(param)`.

## Optional keyword arguments
- `dealias`: dealiasing with `1/3` Orlicz rule if `true` or no dealiasing if `false` (by default);
- `label`: a label for future references (default is `"deep quadratic"`);
- `verbose`: prints information if `true` (default is `true`).

# Return values
Generate necessary ingredients for solving an initial-value problem via `solve!`:
1. a function `DeepQuadratic.f!` to be called in explicit time-integration solvers;
2. a function `DeepQuadratic.mapto` which from `(η,v)` of type `InitialData` provides the raw data matrix on which computations are to be executed;
3. a function `DeepQuadratic.mapfro` which from such data matrix returns the Tuple of real vectors `(η,v)`, where
    - `η` is the surface deformation;
    - `v` is given by `∂t η = - ∂x v`.

"""
mutable struct DeepQuadratic <: AbstractModel

	label   :: String
	f!		:: Function
	mapto	:: Function
	mapfro	:: Function

    function DeepQuadratic( param::NamedTuple;
							dealias=false, label="deep quadratic", verbose=true )

		# Set up
		ϵ = param.ϵ

		if verbose   # Print information
			info = "Build the deep quadratic model.\n"
			info *= "Steepness parameter ϵ=$(param.ϵ) (infinite depth case).\n"
			if dealias == true || dealias == 1
				info *= "Dealiasing with Orszag’s 3/2 rule. "
			else
				info *= "No dealiasing. "
			end
			info *= "\nShut me up with keyword argument `verbose = false`."
			@info info
			@warn "You should provide the initial value for v by ∂t η = - ∂x v.\n\
		It equals the derivative of the trace of the velocity potential \
		used for water waves only when they are null."
		end

		# Pre-allocate useful data
		mesh  = Mesh(param)
		x = mesh.x
		k = mesh.k
        Γ = abs.(k)
        Dx    =  1im * k           	# Differentiation
        H     = -1im * sign.(k)    	# Hilbert transform
		if dealias == true || dealias == 1
			Π⅔    = Γ .< (mesh.kmax-mesh.kmin)/3 	# Dealiasing low-pass filter
		else
			Π⅔    = zero(Γ) .+ 1     		# No dealiasing (Π⅔=Id)
		end

        hnew = zeros(Complex{Float64}, mesh.N)
        unew = zeros(Complex{Float64}, mesh.N)

        I₁ = zeros(Complex{Float64}, mesh.N)
        I₂ = zeros(Complex{Float64}, mesh.N)
        I₃ = zeros(Complex{Float64}, mesh.N)

		# Evolution equations are ∂t U = f(U)
		function f!(U)

			hnew .=ifft(U[:,1]);
			I₁ .=H.*fft(ifft(Γ.*U[:,2]).^2);
			I₂ .=fft(hnew.*ifft(Dx.*U[:,1])) ;
			I₃ .=H.*fft(hnew.*ifft(Γ.*U[:,1]));
			hnew .= U[:,1] ;
			U[:,1] .= -(Dx.*U[:,2]) ;
			U[:,2] .= H.*hnew+ϵ*Π⅔.*(I₁-I₂-I₃) ;

		end

		# Build raw data from physical data.
		# Discrete Fourier transform with, possibly, dealiasing.
		# This function is correct only when the initial data for v is zero.
		function mapto(data::InitialData)
			[Π⅔ .* fft(data.η(x)) Π⅔ .*fft(data.v(x))]
		end

		# Return `(η,v)`, where
		# - `η` is the surface deformation;
		# - `v` is detemined by ∂t η + ∂x v = 0.
		# Inverse Fourier transform and takes the real part.
		function mapfro(U)
			real(ifft(U[:,1])),real(ifft(U[:,2]))
		end

		new(label, f!, mapto, mapfro )

    end
end
