export WWn

@doc raw"""
    WWn(param; kwargs...)

Define an object of type `AbstractModel` in view of solving the initial-value problem for
the water waves expansion proposed by [Dommermuth and Yue](@cite DommermuthYue1987),
[West et al.](@cite WestBruecknerJandaEtAl1987), [Craig and Sulem](@cite CraigSulem1993)
(see also the account by [Choi](@cite Choi2019))
with the "rectification" method proposed by [Duchêne and Melinand](@cite DucheneMelinand2024). 

In the case of the quadratic model (`n=2`) the equations read
```math
  \left\{\begin{array}{l}
  ∂_tη-\tfrac{1}{\sqrtμ ν} T^μv  + \tfrac{ϵ}{ν} ∂_x\big((J^δη) v +  T^μ((J^δη) T^μ v)\big) =0,\\[1ex]
  ∂_tv+∂_xη+\frac{ϵ}{2ν}∂_xJ^δ\big( v^2-(T^μv)^2\big)=0,
  \end{array}\right.
```
where ``η`` is the surface deformation, ``v=∂_xψ`` is the derivative of the trace of the velocity potential at the surface, and
```math
T^μ=-{\rm i}\tanh(\sqrtμ D)
```
is the Fourier multiplier sometimes called "Tilbert transform"
(related to the [Hilbert transform](https://en.wikipedia.org/wiki/Hilbert_transform#Relationship_with_the_Fourier_transform), the latter arising in the infinite layer configuration, ``μ=∞``),
and finally ``J^δ=J_0(δD)`` is the rectification operator: typically  ``J_0(k)=\min(1,1/|k|)`` and the parameter ``δ`` can be freely chosen, but is typically of the size of ``\tfrac{ϵ}{ν}``.

# Argument
`param` is of type `NamedTuple` and must contain
- dimensionless parameters `ϵ` (nonlinearity) and `μ` (dispersion);
- optionally, `ν` the shallow/deep water scaling factor. By default, `ν=1` if `μ≦1` and `ν=1/√μ` otherwise. Set the infinite-layer case if `ν=0`, or `μ=Inf`.
- numerical parameters to construct the mesh of collocation points, if `mesh` is not provided as a keyword argument.

## Optional keyword arguments
- `IL`: Set the infinite-layer case if `IL=true` (or `μ=Inf`, or `ν=0`), in which case `ϵ` is the steepness parameter. Default is `false`.
- `n :: Int`: the order of the expansion; linear system if `1`, quadratic if `2`, cubic if `3`, quartic if `4` (default and other values yield `2`);
- `δ` and `m`: parameters of the rectifier operator, set as `k->min(1,|δ*k|^m)` or `k->min(1,|δ*k|^m[1]*exp(1-|δ*k|^m[2]))` if `m` is a couple
(by default is `δ=0`, i.e. no regularization and `m=-1`. Notice `m=-Inf` and `δ>0` yields a cut-off filter);
- `mesh`: the mesh of collocation points. By default, `mesh = Mesh(param)`;
- `ktol`: tolerance of the low-pass Krasny filter (default is `0`, i.e. no filtering);
- `dealias`: dealiasing with Orszag rule `1-dealias/(dealias+2)` (default is `0`, i.e. no dealiasing);
- `label`: a label for future references (default is `"WWn"` with `n` the order of the expansion);

# Return values
Generate necessary ingredients for solving an initial-value problem via `solve!`:
1. a function `WWn.f!` to be called in explicit time-integration solvers (also `WWn.f1!` and `WWn.f2!` for the symplectic Euler solver);
Additionnally, two functions `WaterWaves.f1!` and `WaterWaves.f2!` for symplectic solvers;
2. a function `WWn.mapto` which from `(η,v)` of type `InitialData` provides the raw data matrix on which computations are to be executed;
3. a function `WWn.mapfro` which from such data matrix returns the Tuple of real vectors `(η,v,x)`, where
  - `η` is the values of surface deformation at collocation points `x`;
  - `v` is the derivative of the trace of the velocity potential at `x`.

"""
mutable struct WWn <: AbstractModel

	label   :: String
	f!		:: Function
	f1!		:: Function
	f2!		:: Function
	mapto	:: Function
	mapfro	:: Function
	info 	:: String

    function WWn(param::NamedTuple;
							mesh = Mesh(param),
							IL	    = false,
							n		= 2,
							δ		= 0,
							m		= -1,
							ktol	= 0,
							dealias	= 0,
							label	= nothing
							)

		# Set up
		if !(n in [1,2,3,4]) n=2 end
		if isnothing(label)  label = "WW$n" end
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
		info = "Spectral model of order $n.\n"
		if IL == true
			info *= "├─Steepness parameter ϵ=$ϵ (infinite depth case).\n"
		else
			info *= "├─Shallowness parameter μ=$μ, nonlinearity parameter ϵ=$ϵ, \
					scaling parameter ν=$nu.\n"
		end
		if δ != 0
			info *= "├─Rectifier with strength δ=$δ and order m=$m.\n"
		end
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

		# Pre-allocate useful data
		k = mesh.k
		x = mesh.x

		if length(m) == 1
			if m==-Inf
				rectifier = k -> k<1
			else
				rectifier = k -> min(1,abs(k)^m)
			end
		else
			rectifier = k -> min(1,abs(k)^(m[1])*exp(1-abs(k)^m[2]))
		end
		Jδ = rectifier.(δ*k)   # regularizing rectifier
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
		η = copy(z) ; v = copy(z) ; Lphi = copy(z) ; LzLphi = copy(z) ; dxv = copy(z) ;
		fftη = copy(z) ; fftv = copy(z) ; Q = copy(z) ; R = copy(z) ;

		# Build raw data from physical data.
		# Discrete Fourier transform with, possibly, dealiasing and Krasny filter.
		function mapto(data::InitialData)
			U = [Π⅔ .* fft(data.η(x)), Π⅔ .*fft(data.v(x))]
			for u in U u[ abs.(u).< ktol ].=0 end
			return U
		end

		# Reconstruct physical variables from raw data
		# Return `(η,v,x)`, where
		# - `η` is the surface deformation;
		# - `v` is the derivative of the trace of the velocity potential;
		# - `x` is the vector of collocation points
		function mapfro(U)
			real( ifft(U[1]) ),real( ifft(U[2]) ), mesh.x
		end

		# Evolution equations are ∂t U = f(U)
		function f!(U)

			fftη .= U[1]
		    fftv .= U[2]
			Q .= Tμ.*fftv
			R .= -fftη*ν

			if n >= 2
				η  .= ifft(Jδ.*U[1])
				v  .= ifft(U[2])
				Lphi .= -ifft(Q)
				Q += -ϵ *∂ₓ.*fft(η.*ifft(fftv)) .+ ϵ*G₀.*fft(η.*Lphi)
				R += ϵ/2*Jδ.*fft(-v.^2 .+ Lphi.^2)
			end
			if n >= 3
				LzLphi .= ifft(-G₀.* fft(η.*Lphi))
				dxv .= ifft(∂ₓ.*fftv)
				Q += ϵ^2*G₀.*fft(η.* LzLphi + 1/2 * η.^2 .* dxv ) .- ϵ^2*∂ₓ.*∂ₓ.*fft( 1/2 * η.^2  .*Lphi)
				R += ϵ^2*Jδ.*fft(Lphi .* ( LzLphi .+ η.* dxv ) )
			end
			if n >= 4
				Q += ϵ^3 * G₀.*fft(η.*ifft(-G₀.* fft(η.*LzLphi + 1/2 * η.^2 .* dxv ) )
							.+ 1/2 * η.^2 .* ifft(∂ₓ.*∂ₓ.*fft( η  .* Lphi ) )
							.- 1/6 * η.^3 .* ifft(∂ₓ.*∂ₓ.*fft( Lphi ) ) ) .-
					   ϵ^3 * ∂ₓ.*∂ₓ.*fft( 1/2 * η.^2  .*LzLphi .+ 1/3 * η.^3  .* dxv )
				R += ϵ^3 * Jδ.*fft( Lphi .*  ifft(-G₀.* fft(η.*LzLphi + 1/2 * η.^2 .* dxv ) )
						.+ 1/2* (LzLphi .+ η .* dxv ).^2
						.+ 1/2* η.* Lphi.^2 .* ifft(∂ₓ.*∂ₓ.* fftη)
						.- 1/2* (η.^2).* (ifft(∂ₓ .* fft(Lphi))).^2 ) .+
						1/4*ϵ^3 * ∂ₓ.*∂ₓ.*fft((η .* Lphi).^2)
			end
		   	U[1] .= Π⅔.*Q/sqrt(μ)/ν
		   	U[2] .= Π⅔.*∂ₓ.*R/sqrt(μ)/ν
			for u in U u[ abs.(u).< ktol ].=0 end

		end

		# Evolution equations are ∂t (U1,U2) = (f1(U1,U2) , f2(U1,U2))
		function f1!(U1,U2)

			fftη .= U1
		    fftv .= U2
			Q .= Tμ.*fftv

			if n >= 2
				η  .= ifft(Jδ.*U1)
				v  .= ifft(U2)
				Lphi .= -ifft(Q)
				Q += -ϵ *∂ₓ.*fft(η.*ifft(fftv)) .+ ϵ*G₀.*fft(η.*Lphi)
			end
			if n >= 3
				LzLphi .= ifft(-G₀.* fft(η.*Lphi))
				dxv .= ifft(∂ₓ.*fftv)
				Q += ϵ^2*G₀.*fft(η.* LzLphi + 1/2 * η.^2 .* dxv ) .-
						ϵ^2*∂ₓ.*∂ₓ.*fft( 1/2 * η.^2  .*Lphi)
			end
			if n >= 4
				Q += ϵ^3 * G₀.*fft(η.*ifft(-G₀.* fft(η.*LzLphi + 1/2 * η.^2 .* dxv ) )
							.+ 1/2 * η.^2 .* ifft(∂ₓ.*∂ₓ.*fft( η  .* Lphi ) )
							.- 1/6 * η.^3 .* ifft(∂ₓ.*∂ₓ.*fft( Lphi ) ) ) .-
					   ϵ^3 * ∂ₓ.*∂ₓ.*fft( 1/2 * η.^2  .*LzLphi .+ 1/3 * η.^3  .* dxv )
			end
		   	U1 .= Π⅔.*Q/sqrt(μ)/ν
			U1[abs.(U1).< ktol ].=0

		end

		function f2!(U1,U2)

			fftη .= U1
		    fftv .= U2
			Q .= Tμ.*fftv
			R .= -fftη*ν

			if n >= 2
				η  .= ifft(Jδ.*U1)
				v  .= ifft(U2)
				Lphi .= -ifft(Q)
				R += ϵ/2*Jδ.*fft(-v.^2 .+ Lphi.^2)
			end
			if n >= 3
				LzLphi .= ifft(-G₀.* fft(η.*Lphi))
				dxv .= ifft(∂ₓ.*fftv)
				R += ϵ^2*Jδ.*fft(Lphi .* ( LzLphi .+ η.* dxv ) )
			end
			if n >= 4
				R += ϵ^3 * Jδ.*fft( Lphi .*  ifft(-G₀.* fft(η.*LzLphi + 1/2 * η.^2 .* dxv ) )
						.+ 1/2* (LzLphi .+ η .* dxv ).^2
						.+ 1/2* η.* Lphi.^2 .* ifft(∂ₓ.*∂ₓ.* fftη)
						.- 1/2* (η.^2).* (ifft(∂ₓ .* fft(Lphi))).^2 ) .+
						1/4*ϵ^3 * ∂ₓ.*∂ₓ.*fft((η .* Lphi).^2)
			end
		   	U2 .= Π⅔.*∂ₓ.*R/sqrt(μ)/ν
			U2[ abs.(U2).< ktol ].=0

		end

        new(label, f!, f1!, f2! , mapto, mapfro, info )
    end
end
