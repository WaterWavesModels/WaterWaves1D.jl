export PseudoSpectral,mapto,mapfro

"""
    PseudoSpectral(params)
	Uses a PseudoSpectral expansion (WWn), see West et al., Craig-Sulem, etc.
	The operators are damped in a way that preserves the asymptotic.
	params must contain "order", an integer in {1,2,3}. See Choi in RIMS Workshop for explicit expression order up to 5

"""
mutable struct PseudoSpectral <: AbstractModel

    label   :: String
	datasize:: Int
	μ 		:: Float64
	ϵ 		:: Float64
	n 		:: Int
	x   	:: Array{Float64,1}
    F₀   	:: Array{Float64,1}
	G₀   	:: Array{Float64,1}
    ∂ₓ      :: Array{Complex{Float64},1}
    Π⅔      :: BitArray{1}
	η    	:: Vector{Float64}
    fftη    :: Vector{Complex{Float64}}
    fftv    :: Vector{Complex{Float64}}
	Q	    :: Vector{Complex{Float64}}
    R	    :: Vector{Complex{Float64}}

    function PseudoSpectral(param::NamedTuple)

		label = "Pseudo Spectral"
		datasize = 2
		μ 	= param.μ
		ϵ 	= param.ϵ
		n 	= param.order
		mesh = Mesh(param)
		x = mesh.x
		k = copy(mesh.k)
		cutoff = k -> tanh.(-abs.(k).+1 ./abs.(k))
		Π = cutoff(10*ϵ*k)   # regularisation cutoff
		F₀ 	= tanh.(sqrt(μ)*abs.(mesh.k))./(sqrt(μ)*abs.(mesh.k))
		F₀[1] 	= 1
		G₀ 	= 1/(sqrt(μ))*Π.*abs.(mesh.k).*tanh.(sqrt(μ)*abs.(mesh.k))
		∂ₓ	=  1im * Π.* mesh.k            # Differentiation
        Π⅔ 	= abs.(mesh.k) .< mesh.kmax * 2/3 # Dealiasing low-pass filter
		η = zeros(Float64, mesh.N)
		z = zeros(Complex{Float64}, mesh.N)
		fftη = copy(z) ; fftv = copy(z) ; Q = copy(z) ; R = copy(z) ;

        new(label, datasize, μ, ϵ, n, x, F₀, G₀, ∂ₓ, Π⅔, η, fftη, fftv, Q, R )
    end
end


function (m::PseudoSpectral)(U::Array{Complex{Float64},2})


	m.fftη .= U[:,1]
    m.fftv .= U[:,2]
	m.η  .= real(ifft(U[:,1]))
   	m.Q .= -m.∂ₓ.*m.F₀.*m.fftv #-Lϕ=G₀ϕ
	m.R .= -m.fftη
	# attention,  m.G₀=-L dans Choi
	if m.n >= 2
		m.Q += -m.ϵ *m.∂ₓ.*fft(m.η.*ifft(m.fftv)) .- m.ϵ*m.G₀.*fft(m.η.*ifft(m.Q))
		m.R += m.ϵ/2*fft(-ifft(m.fftv).^2 .+ ifft(m.Q).^2)
	end
	if m.n >= 3
		m.Q += m.ϵ^2*m.G₀.*fft(m.η.*ifft(m.G₀.* fft(m.η.*ifft(m.Q)) ) + 1/2 * m.η.^2 .* ifft(m.∂ₓ.*m.fftv) ) .+ m.ϵ^2*m.∂ₓ.*m.∂ₓ.*fft( 1/2 * m.η.^2  .*ifft(m.Q))
		m.R += m.ϵ^2*fft(ifft(-m.Q) .*  ifft(m.G₀.*fft(m.η.*ifft(m.Q)) .+ m.η.* ifft(m.∂ₓ.*m.fftv)) )
	end
   	U[:,1] .= m.Q
   	U[:,2] .= m.∂ₓ.*m.R

end

"""
    mapto(PseudoSpectral, data)

"""
function mapto(m::PseudoSpectral, data::InitialData)

	[m.Π⅔ .* fft(data.η(m.x)) m.Π⅔ .*fft(data.v(m.x))]

end

"""
    mapfro(PseudoSpectral, data)

"""
function mapfro(m::PseudoSpectral,
	       datum::Array{Complex{Float64},2})
		   real(ifft(datum[:,1])),real(ifft(datum[:,2]))
end
