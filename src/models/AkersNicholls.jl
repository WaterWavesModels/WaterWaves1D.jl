export AkersNicholls_fast,AkersNicholls

"""
    AkersNicholls_fast(param; kwargs...)

Same as [`AkersNicholls`](@ref), but faster.
"""
mutable struct AkersNicholls_fast <: AbstractModel

	label   :: String
	f!		:: Function
	mapto	:: Function
	mapfro	:: Function
	info    :: String

    function AkersNicholls_fast( param::NamedTuple;
									mesh = Mesh(param),
									IL	    = false,
									dealias=false,
									label="Akers-Nicholls" )

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

		# Print information
		info = "Akers-Nicholls model.\n"
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
			Lμ 	= abs.(k);Lμ[1]=1;
		else
			Tμ 	= -1im* sign.(k) .* tanh.(sqrt(μ)*abs.(k))
			Lμ 	= sqrt(μ)*abs.(k)./tanh.(sqrt(μ)*abs.(k));Lμ[1]=1;
		end
		∂ₓ	=  1im * sqrt(μ)* k            # Differentiation
		Tμ∂ₓ= Tμ .* ∂ₓ
		if dealias == true || dealias == 1
			Π⅔    = abs.(k) .< (mesh.kmax-mesh.kmin)/3 	# Dealiasing low-pass filter
		else
			Π⅔    = zero(k) .+ 1     		# No dealiasing (Π⅔=Id)
		end

        ζ = zeros(Complex{Float64}, mesh.N)
        unew = zeros(Complex{Float64}, mesh.N)

        I₁ = zeros(Complex{Float64}, mesh.N)
        I₂ = zeros(Complex{Float64}, mesh.N)
        I₃ = zeros(Complex{Float64}, mesh.N)

        Px  = plan_fft(ζ; flags = FFTW.MEASURE)

		# Evolution equations are ∂t U = f(U)
		function f!(U)

		    ldiv!(ζ, Px , U[1])

		    I₁  .= U[2] .* Lμ
		    ldiv!(unew, Px , I₁)
		    unew .^= 2
		    mul!(I₁, Px , unew)
		    I₁ .*= Tμ

		    I₂  .= U[1] .* ∂ₓ
		    ldiv!(unew, Px , I₂)
		    unew .*= ζ
		    mul!(I₂, Px , unew)

		    I₃  .= U[1] .* Tμ∂ₓ
		    ldiv!(unew, Px, I₃)
		    unew .*= ζ
		    mul!(I₃ , Px , unew)
		    I₃ .*= Tμ

		    ζ  .= .- U[2] .* ∂ₓ

		    I₁ .*= ν
			I₁ .-= I₂
		    I₁ .-= I₃
		    I₁ .*= Π⅔
		    I₁ .*= ϵ

		    U[2] .= (U[1] .* Tμ .+ I₁)/sqrt(μ)/ν
		    U[1] .= ζ/sqrt(μ)

		end

		# Build raw data from physical data.
		function mapto(data::InitialData)
			ζ.=data.η(x);v=data.v(x);
			fftv=fft(v);
			[Π⅔ .* fft(ζ) , Π⅔ .* (fftv./Lμ+ϵ*fft(ζ.*v)+ ϵ*Tμ.*fft(ζ.*ifft(Tμ.*fftv)))/ν]
		end

		# Reconstruct physical variables from raw data
		# Return `(η,v,x)`, where
		# - `η` is the surface deformation;
		# - `v` is the derivative of the trace of the velocity potential;
		# - `x` is the vector of collocation points
		function mapfro(U;n=10)
			ζ.=ifft(U[1]);I₁.=ν*U[2];I₂.=Lμ.*I₁;
			for j=1:n
				I₂.=Lμ.*(I₁-ϵ*Π⅔ .* ( fft(ζ.*ifft(I₂)) + Tμ.*fft(ζ.*ifft(Tμ.*I₂))) )
			end
			real(ζ),real(ifft(I₂)),mesh.x
		end

		new(label, f!, mapto, mapfro, info )

    end
end

@doc raw"""
    AkersNicholls(param; kwargs...)

Define an object of type `AbstractModel` in view of solving the initial-value problem for
the quadratic deep-water model proposed by [AkersNicholls2010](@citet)
and [ChengGranero-BelinchonShkollerEtAl2019](@citet)
```math
  \left\{\begin{array}{l}
  ∂_tη+∂_x m=0,\\[1ex]
  ∂_tm-\tfrac{1}{\sqrtμ ν} T^μ\big(η+\frac{ϵ}{ν}(L^μ m)^2\big)+\frac{ϵ}{ν}\big(η∂_xη+T^μ(η ∂_x T^μ η)\big)=0,
  \end{array}\right.
```
where ``η`` is the surface deformation, ``m=-\frac1{\sqrtμ ν} T^μψ  + \frac{ϵ}{ν} \big(η ∂_xψ +  T^μ(η T^μ ∂_xψ)\big)`` represents the vertically integrated horizontal momentum, and
``T^μ=-{\rm i}\tanh(\sqrtμ D)`` and ``L^μ=\frac{ν\sqrtμ D}{\tanh(\sqrtμ D)}`` are Fourier multipliers.

# Argument
`param` is of type `NamedTuple` and must contain
- dimensionless parameters `ϵ` (nonlinearity) and `μ` (dispersion);
- optionally, `ν` the shallow/deep water scaling factor. By default, `ν=1` if `μ≦1` and `ν=1/√μ` otherwise. Set the infinite-layer case if `ν=0`, or `μ=Inf`.
- numerical parameters to construct the mesh of collocation points, if `mesh` is not provided as a keyword argument.

## Optional keyword arguments
- `mesh`: the mesh of collocation points. By default, `mesh = Mesh(param)`;
- `IL`: Set the infinite-layer case if `IL=true` (or `μ=Inf`, or `ν=0`), in which case `ϵ` is the steepness parameter. Default is `false`.
- `dealias`: dealiasing with `1/3` Orszag rule if `true` or no dealiasing if `false` (by default);
- `label`: a label for future references (default is `"deep quadratic"`);

# Return values
Generate necessary ingredients for solving an initial-value problem via `solve!`:
1. a function `AkersNicholls.f!` to be called in explicit time-integration solvers;
2. a function `AkersNicholls.mapto` which from `(η,v)` of type `InitialData` provides the raw data matrix on which computations are to be executed;
3. a function `AkersNicholls.mapfro` which from such data matrix returns the Tuple of real vectors `(η,v,x)`, where
    - `η` is the values of surface deformation at collocation points `x`;
    - `v` is the derivative of the trace of the velocity potential at `x`.

Consider also [`AkersNicholls_fast`](@ref).
"""
mutable struct AkersNicholls <: AbstractModel

	label   :: String
	f!		:: Function
	mapto	:: Function
	mapfro	:: Function
	info 	:: String

    function AkersNicholls( param::NamedTuple;
							mesh = Mesh(param),
							IL	    = false,
							dealias = false,
							label="Akers-Nicholls")

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

		# Print information
		info = "Akers-Nicholls model.\n"
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
			Lμ 	= abs.(k);Lμ[1]=1;
		else
			Tμ 	= -1im* sign.(k) .* tanh.(sqrt(μ)*abs.(k))
			Lμ 	= sqrt(μ)*abs.(k)./tanh.(sqrt(μ)*abs.(k));Lμ[1]=1;
		end
		∂ₓ	=  1im * sqrt(μ)* k            # Differentiation
		Tμ∂ₓ= Tμ .* ∂ₓ

		if dealias == true || dealias == 1
			Π⅔    = abs.(k) .< (mesh.kmax-mesh.kmin)/3 	# Dealiasing low-pass filter
		else
			Π⅔    = zero(k) .+ 1     		# No dealiasing (Π⅔=Id)
		end

        ζ = zeros(Complex{Float64}, mesh.N)
        unew = zeros(Complex{Float64}, mesh.N)

        I₁ = zeros(Complex{Float64}, mesh.N)
        I₂ = zeros(Complex{Float64}, mesh.N)
        I₃ = zeros(Complex{Float64}, mesh.N)

		# Evolution equations are ∂t U = f(U)
		function f!(U)

			ζ .=ifft(U[1]);
			I₁ .=Tμ.*fft(ifft(Lμ.*U[2]).^2);
			# Tricomi identity : the above equals the two commented lines below
			# I₁ .=Tμ.*fft(ifft(Lμ.*U[2]).^2 .+ ifft(∂ₓ.*U[2]).^2)/2;
			# I₁ .-=fft(ifft(Lμ.*U[2]).*ifft(∂ₓ.*U[2]));

			I₂ .=fft(ζ.*ifft(∂ₓ.*U[1])) ;
			I₃ .=Tμ.*fft(ζ.*ifft(Tμ∂ₓ.*U[1]));
			ζ .= U[1] ;
			U[1] .= -∂ₓ.*U[2]/sqrt(μ) ;
			U[2] .= (Tμ.*ζ+ϵ*Π⅔.*(ν*I₁-I₂-I₃))/sqrt(μ)/ν ;

		end

		# Build raw data from physical data.
		function mapto(data::InitialData)
			ζ.=data.η(x);v=data.v(x);
			fftv=fft(v);
			[Π⅔ .* fft(ζ),  Π⅔ .* (fftv./Lμ+ϵ*fft(ζ.*v)+ ϵ*Tμ.*fft(ζ.*ifft(Tμ.*fftv)))/ν]
		end

		# Reconstruct physical variables from raw data
		# Return `(η,v,x)`, where
		# - `η` is the surface deformation;
		# - `v` is the derivative of the trace of the velocity potential;
		# - `x` is the vector of collocation points
		function mapfro(U;n=10)
			ζ.=ifft(U[1]);I₁.=ν*U[2];I₂.=Lμ.*I₁;
			for j=1:n
				I₂.=Lμ.*(I₁-ϵ*Π⅔ .* ( fft(ζ.*ifft(I₂)) + Tμ.*fft(ζ.*ifft(Tμ.*I₂))) )
			end
			real(ζ),real(ifft(I₂)),mesh.x
		end

		new(label, f!, mapto, mapfro, info )

    end
end
