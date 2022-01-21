export WWn

"""
    WWn(param;kwargs)

Define an object of type `AbstractModel` in view of solving the initial-value problem for
the water waves expansion proposed by Dommermuth and Yue (doi:10.1017/s002211208700288x), West et al. (doi:10.1029/jc092ic11p11803), Craig and Sulem (doi:10.1006/jcph.1993.1164)
(see also the account by Choi : Fifth-order nonlinear spectral model for surface gravity waves: From pseudo-spectral to spectral formulations, RIMS Kokyuroku, 2019)
with the "rectification" method proposed by Duchêne and Melinand.

# Argument
`param` is of type `NamedTuple` and must contain
- dimensionless parameters `ϵ` (nonlinearity) and `μ` (dispersion);
- numerical parameters to construct the mesh of collocation points as `mesh = Mesh(param)`

## Optional keyword arguments
- `ν`: shallow/deep water multiplication factor. By default, `ν=1` if `μ≦1` and `ν=1/√μ` otherwise. Set the infinite-layer case if `ν=0` (or `μ=Inf`).
- `n :: Int`: the order of the expansion; linear system if `1`, quadratic if `2`, cubic if `3`, quartic if `4` (default and other values yield `2`);
- `δ` and `m`: parameters of the rectifier operator, set as `k->min(1,|δ*k|^m)` or `k->min(1,|δ*k|^m[1]*exp(1-|δ*k|^m[2]))` if `m` is a couple
(by default is `δ=0`, i.e. no regularization and `m=-1`. Notice `m=-Inf` and `δ>0` yields a cut-off filter);
- `ktol`: tolerance of the low-pass Krasny filter (default is `0`, i.e. no filtering);
- `dealias`: dealiasing with Orlicz rule `1-dealias/(dealias+2)` (default is `0`, i.e. no dealiasing);
- `label`: a label for future references (default is `"WWn"` with `n` the order of the expansion);
- `verbose`: prints information if `true` (default is `true`).

# Return values
Generate necessary ingredients for solving an initial-value problem via `solve!`:
1. a function `WWn.f!` to be called in explicit time-integration solvers (also `WWn.f1!` and `WWn.f2!` for the symplectic Euler solver);
2. a function `WWn.mapto` which from `(η,v)` of type `InitialData` provides the raw data matrix on which computations are to be executed;
3. a function `WWn.mapfro` which from such data matrix returns the Tuple of real vectors `(η,v)`, where
    - `η` is the surface deformation;
    - `v` is the derivative of the trace of the velocity potential.

"""
mutable struct WWn <: AbstractModel

	label   :: String
	f!		:: Function
	f1!		:: Function
	f2!		:: Function
	mapto	:: Function
	mapfro	:: Function

    function WWn(param::NamedTuple;
							ν		= nothing,
							n		= 2,
							δ		= 0,
							m		= -1,
							ktol	= 0,
							dealias	= 0,
							label	= nothing,
							verbose	= true)

		if !(n in [1,2,3,4]) n=2 end
		if label == nothing  label = "WW$n" end
		if verbose @info "Build the spectral model with nonlinearity of order $n" end

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

		if length(m) == 1
			rectifier = k -> min(1,abs(k)^m)
		else
			rectifier = k -> min(1,abs(k)^(m[1])*exp(1-abs(k)^m[2]))
		end
		Π = rectifier.(δ*k)   # regularizing rectifier
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
		η = copy(z) ; v = copy(z) ; Lphi = copy(z) ; LzLphi = copy(z) ; dxv = copy(z) ;
		fftη = copy(z) ; fftv = copy(z) ; Q = copy(z) ; R = copy(z) ;

		# Discrete Fourier transform with, possibly, dealiasing and Krasny filter.
		function mapto(data::InitialData)
			U = [Π⅔ .* fft(data.η(x)) Π⅔ .*fft(data.v(x))]
			U[ abs.(U).< ktol ].=0
			return U
		end

		# Return `(η,v)`, where
		# - `η` is the surface deformation;
		# - `v` is the derivative of the trace of the velocity potential.
		# Inverse Fourier transform and takes the real part.
		function mapfro(U)
			real( ifft(U[:,1]) ),real( ifft(U[:,2]) )
		end

		# Evolution equations are ∂t U = f!(U)
		function f!(U)

			fftη .= U[:,1]
		    fftv .= U[:,2]
			Q .= -∂ₓF₀.*fftv
			R .= -fftη*ν

			if n >= 2
				η  .= ifft(Π.*U[:,1])
				v  .= ifft(U[:,2])
				Lphi .= -ifft(Q)
				Q += -ϵ *∂ₓ.*fft(η.*ifft(fftv)) .+ ϵ*G₀.*fft(η.*Lphi)
				R += ϵ/2*Π.*fft(-v.^2 .+ Lphi.^2)
			end
			if n >= 3
				LzLphi .= ifft(-G₀.* fft(η.*Lphi))
				dxv .= ifft(∂ₓ.*fftv)
				Q += ϵ^2*G₀.*fft(η.* LzLphi + 1/2 * η.^2 .* dxv ) .- ϵ^2*∂ₓ.*∂ₓ.*fft( 1/2 * η.^2  .*Lphi)
				R += ϵ^2*Π.*fft(Lphi .* ( LzLphi .+ η.* dxv ) )
			end
			if n >= 4
				Q += ϵ^3 * G₀.*fft(η.*ifft(-G₀.* fft(η.*LzLphi + 1/2 * η.^2 .* dxv ) )
							.+ 1/2 * η.^2 .* ifft(∂ₓ.*∂ₓ.*fft( η  .* Lphi ) )
							.- 1/6 * η.^3 .* ifft(∂ₓ.*∂ₓ.*fft( Lphi ) ) ) .-
					   ϵ^3 * ∂ₓ.*∂ₓ.*fft( 1/2 * η.^2  .*LzLphi .+ 1/3 * η.^3  .* dxv )
				R += ϵ^3 * Π.*fft( Lphi .*  ifft(-G₀.* fft(η.*LzLphi + 1/2 * η.^2 .* dxv ) )
						.+ 1/2* (LzLphi .+ η .* dxv ).^2
						.+ 1/2* η.* Lphi.^2 .* ifft(∂ₓ.*∂ₓ.* fftη)
						.- 1/2* (η.^2).* (ifft(∂ₓ .* fft(Lphi))).^2 ) .+
						1/4*ϵ^3 * ∂ₓ.*∂ₓ.*fft((η .* Lphi).^2)
			end
		   	U[:,1] .= Π⅔.*Q/sqrt(μ)/ν
		   	U[:,2] .= Π⅔.*∂ₓ.*R/sqrt(μ)/ν
			U[abs.(U).< ktol ].=0

		end

		# Evolution equations are ∂t U = (f1!(U) f2!(U))
		function f1!(U1,U2)

			fftη .= U1
		    fftv .= U2
			Q .= -∂ₓF₀.*fftv

			# attention,  G₀=-L dans Choi
			if n >= 2
				η  .= ifft(Π.*U1)
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
			Q .= -∂ₓF₀.*fftv
			R .= -fftη*ν

			# attention,  G₀=-L dans Choi
			if n >= 2
				η  .= ifft(Π.*U1)
				v  .= ifft(U2)
				Lphi .= -ifft(Q)
				R += ϵ/2*Π.*fft(-v.^2 .+ Lphi.^2)
			end
			if n >= 3
				LzLphi .= ifft(-G₀.* fft(η.*Lphi))
				dxv .= ifft(∂ₓ.*fftv)
				R += ϵ^2*Π.*fft(Lphi .* ( LzLphi .+ η.* dxv ) )
			end
			if n >= 4
				R += ϵ^3 * Π.*fft( Lphi .*  ifft(-G₀.* fft(η.*LzLphi + 1/2 * η.^2 .* dxv ) )
						.+ 1/2* (LzLphi .+ η .* dxv ).^2
						.+ 1/2* η.* Lphi.^2 .* ifft(∂ₓ.*∂ₓ.* fftη)
						.- 1/2* (η.^2).* (ifft(∂ₓ .* fft(Lphi))).^2 ) .+
						1/4*ϵ^3 * ∂ₓ.*∂ₓ.*fft((η .* Lphi).^2)
			end
		   	U2 .= Π⅔.*∂ₓ.*R/sqrt(μ)/ν
			U2[abs.(U2).< ktol ].=0

		end



        new(label, f!, f1!, f2! , mapto, mapfro )
    end
end
