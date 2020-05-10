export PseudoSpectral
using FFTW
"""
    PseudoSpectral(param;kwargs)

Defines an object of type `AbstractModel` in view of solving the initial-value problem for
the modified water waves expansion proposed by West et al., Craig-Sulem, etc.

# Argument
`param` is of type `NamedTuple` and must contain
- dimensionless parameters `ϵ` (nonlinearity) and `μ` (dispersion);
- numerical parameters to construct the mesh of collocation points as `mesh = Mesh(param)`
E.g. `param = ( μ  = 0.1, ϵ  = 1, N  = 2^9, L  = 10*π)`

## Keywords
- `order :: Int`: the order of the expansion; linear system if `1`, quadratic if `2`, cubic if `3`, quartic if `4` (default and other values yield `2`);
- `lowpass :: Real`: applies a low-pass filter on frequencies above `lowpass` (default is `0`, i.e. no low-pass filter);
- `dealias :: Int`: dealiasing with Orlicz rule `1-dealias/(dealias+2)` (default is `0`, i.e. no dealiasing).

# Return values
This generates
1. a function `PseudoSpectral` to be called in the time-integration solver;
2. a function `mapto` which from `(η,v)` of type `InitialData` provides the raw data matrix on which computations are to be executed.
3. a function `mapfro` which from such raw data matrix returns the Tuple of real vectors `(η,v)`, where

    - `η` is the surface deformation;
    - `v` is a tangential surface velocity (the derivative of the trace of the velocity potential at the surface).
"""
mutable struct PseudoSpectral <: AbstractModel

    label   :: String
	mapto	:: Function
	mapfro	:: Function
	f!		:: Function
	μ 		:: Float64
	ϵ 		:: Float64
	n 		:: Int
	x   	:: Vector{Float64}
    F₀   	:: Vector{Float64}
	G₀   	:: Vector{Float64}
    ∂ₓ      :: Vector{Complex{Float64}}
	Π   	:: Vector{Float64}
    Π⅔      :: BitArray{1}
	η	    :: Vector{Complex{Float64}}
    v	    :: Vector{Complex{Float64}}
	fftη    :: Vector{Complex{Float64}}
    fftv    :: Vector{Complex{Float64}}
	Q	    :: Vector{Complex{Float64}}
    R	    :: Vector{Complex{Float64}}
	Lphi	:: Vector{Complex{Float64}}
	LzLphi	:: Vector{Complex{Float64}}
	dxv		:: Vector{Complex{Float64}}


    function PseudoSpectral(param::NamedTuple;order=2::Int,lowpass=0::Real,dealias=0::Int,)
		if order in [1,2,3,4]
			@info string("solving system with nonlinearity of order ", order)
		else
			order = 2
			@info string("solving system with nonlinearity of order ", order)
		end
		label = string("WW",order)
		datasize = 2
		μ 	= param.μ
		ϵ 	= param.ϵ
		n 	= order
		mesh = Mesh(param)
		x = mesh.x
		k = copy(mesh.k)
		cutoff = k -> (1 + tanh(-abs(k)+1/abs(k)))/2
		Π = cutoff.(lowpass*k)   # regularisation cutoff
		K = mesh.kmax * (1-dealias/(2+dealias))
		Π⅔ 	= abs.(mesh.k) .<= K # Dealiasing low-pass filter
		if dealias == 0
			@info "no dealiasing"
			Π⅔ 	= ones(size(mesh.k))
		else
			@info "dealiasing"
		end
		F₀ 	= tanh.(sqrt(μ)*abs.(mesh.k))./(sqrt(μ)*abs.(mesh.k))
		F₀[1] 	= 1
		G₀ 	= sqrt(μ)*abs.(mesh.k).*tanh.(sqrt(μ)*abs.(mesh.k))
		∂ₓ	=  1im * sqrt(μ)* mesh.k            # Differentiation
		z = zeros(Complex{Float64}, mesh.N)
		η = copy(z) ; v = copy(z) ; Lphi = copy(z) ; LzLphi = copy(z) ; dxv = copy(z) ;
		fftη = copy(z) ; fftv = copy(z) ; Q = copy(z) ; R = copy(z) ;

		function mapto(data::InitialData)

			[Π⅔ .* fft(data.η(x)) Π⅔ .*fft(data.v(x))]

		end

		function mapfro(U)

				   real(ifft(U[:,1])),real(ifft(U[:,2]))
		end

		function f!(U) 		# utile uniquement pour solve2!

			fftη .= U[:,1]
		    fftv .= U[:,2]
			Q .= -∂ₓ.*F₀.*fftv
			R .= -fftη

			# attention,  G₀=-L dans Choi
			if n >= 2
				η  .= ifft(Π.*U[:,1])
				v  .= ifft(U[:,2])
				Lphi .= -ifft(Q)
				Q += -ϵ *∂ₓ.*fft(η.*ifft(fftv)) .+ ϵ*G₀.*fft(η.*Lphi)
				R += ϵ/2*fft(-v.^2 .+ Lphi.^2)
			end
			if n >= 3
				LzLphi .= ifft(-G₀.* fft(η.*Lphi))
				dxv .= ifft(∂ₓ.*fftv)
				Q += ϵ^2*G₀.*fft(η.* LzLphi + 1/2 * η.^2 .* dxv ) .- ϵ^2*∂ₓ.*∂ₓ.*fft( 1/2 * η.^2  .*Lphi)
				R += ϵ^2*fft(Lphi .* ( LzLphi .+ η.* dxv ) )
			end
			if n >= 4
				Q += ϵ^3 * G₀.*fft(η.*ifft(-G₀.* fft(η.*LzLphi + 1/2 * η.^2 .* dxv ) )
							.+ 1/2 * η.^2 .* ifft(∂ₓ.*∂ₓ.*fft( η  .* Lphi ) )
							.- 1/6 * η.^3 .* ifft(∂ₓ.*∂ₓ.*fft( Lphi ) ) ) .-
					   ϵ^3 * ∂ₓ.*∂ₓ.*fft( 1/2 * η.^2  .*LzLphi .+ 1/3 * η.^3  .* dxv )
				R += ϵ^3 * fft( Lphi .*  ifft(-G₀.* fft(η.*LzLphi + 1/2 * η.^2 .* dxv ) )
						.+ 1/2* (LzLphi .+ η .* dxv ).^2
						.+ 1/2* η.* Lphi.^2 .* ifft(∂ₓ.*∂ₓ.* fftη)
						.- 1/2* (η.^2).* (ifft(∂ₓ .* fft(Lphi))).^2 ) .+
						1/4*ϵ^3 * ∂ₓ.*∂ₓ.*fft((η .* Lphi).^2)
			end
		   	U[:,1] .= Π⅔.*Q/sqrt(μ)
		   	U[:,2] .= Π⅔.*∂ₓ.*Π.*R/sqrt(μ)

		end


        new(label, mapto, mapfro, f!, μ, ϵ, n, x, F₀, G₀, ∂ₓ, Π, Π⅔, η, v, fftη, fftv, Q, R, Lphi, LzLphi, dxv )
    end
end

function (m::PseudoSpectral)(U)  # utile uniquement pour solve!

	m.fftη .= U[:,1]
    m.fftv .= U[:,2]
	m.Q .= -m.∂ₓ.*m.F₀.*m.fftv
	m.R .= -m.fftη

	# attention,  m.G₀=-L dans Choi
	if m.n >= 2
		m.η  .= ifft(m.Π.*U[:,1])
		m.v  .= ifft(U[:,2])
		m.Lphi .= -ifft(m.Q)
		m.Q += -m.ϵ *m.∂ₓ.*fft(m.η.*ifft(m.fftv)) .+ m.ϵ*m.G₀.*fft(m.η.*m.Lphi)
		m.R += m.ϵ/2*fft(-m.v.^2 .+ m.Lphi.^2)
	end
	if m.n >= 3
		m.LzLphi .= ifft(-m.G₀.* fft(m.η.*m.Lphi))
		m.dxv .= ifft(m.∂ₓ.*m.fftv)
		m.Q += m.ϵ^2*m.G₀.*fft(m.η.* m.LzLphi + 1/2 * m.η.^2 .* m.dxv ) .- m.ϵ^2*m.∂ₓ.*m.∂ₓ.*fft( 1/2 * m.η.^2  .*m.Lphi)
		m.R += m.ϵ^2*fft(m.Lphi .* ( m.LzLphi .+ m.η.* m.dxv ) )
	end
	if m.n >= 4
		m.Q += m.ϵ^3 * m.G₀.*fft(m.η.*ifft(-m.G₀.* fft(m.η.*m.LzLphi + 1/2 * m.η.^2 .* m.dxv ) )
					.+ 1/2 * m.η.^2 .* ifft(m.∂ₓ.*m.∂ₓ.*fft( m.η  .* m.Lphi ) )
					.- 1/6 * m.η.^3 .* ifft(m.∂ₓ.*m.∂ₓ.*fft( m.Lphi ) ) ) .-
			   m.ϵ^3 * m.∂ₓ.*m.∂ₓ.*fft( 1/2 * m.η.^2  .*m.LzLphi .+ 1/3 * m.η.^3  .* m.dxv )
		m.R += m.ϵ^3 * fft( m.Lphi .*  ifft(-m.G₀.* fft(m.η.*m.LzLphi + 1/2 * m.η.^2 .* m.dxv ) )
				.+ 1/2* (m.LzLphi .+ m.η .* m.dxv ).^2
				.+ 1/2* m.η.* m.Lphi.^2 .* ifft(m.∂ₓ.*m.∂ₓ.* m.fftη)
				.- 1/2* (m.η.^2).* (ifft(m.∂ₓ .* fft(m.Lphi))).^2 ) .+
				1/4*m.ϵ^3 * m.∂ₓ.*m.∂ₓ.*fft((m.η .* m.Lphi).^2)
	end
   	U[:,1] .= m.Π⅔.*m.Q/sqrt(m.μ)
   	U[:,2] .= m.Π⅔.*m.∂ₓ.*m.Π.*m.R/sqrt(m.μ)

end
