export WhithamBoussinesq,mapto,mapfro

"""
    WhithamBoussinesq(params)
	params must contain a parameter α.
	If α = 1, then the model has been introduced and studied by E. Dinvay and collaborators.
	If α = 1/2, then the model is a quasilinear version.
	If α < 1/2, then expect instabilities stemming from ill-posedness of the model.
"""
mutable struct WhithamBoussinesq <: AbstractModel

    label   :: String
	datasize:: Int
	μ 		:: Float64
	ϵ 		:: Float64
	x   	:: Vector{Float64}
    F₁   	:: Vector{Float64}
	F₂   	:: Vector{Float64}
    ∂ₓ      :: Vector{Complex{Float64}}
    Π⅔      :: BitArray{1}
	η    	:: Vector{Float64}
	v    	:: Vector{Float64}
    fftη    :: Vector{Complex{Float64}}
    fftv    :: Vector{Complex{Float64}}

    function WhithamBoussinesq(param::NamedTuple)

		label = string("Whitham-Boussinesq")
		datasize = 2
		μ 	= param.μ
		ϵ 	= param.ϵ
		mesh = Mesh(param)
		x 	= mesh.x
        F₁ 	= tanh.(sqrt(μ)*abs.(mesh.k))./(sqrt(μ)*abs.(mesh.k))
		F₁[1] 	= 1
		F₂ = F₁.^(param.α)
    	∂ₓ	=  1im * mesh.k            # Differentiation
        Π⅔ 	= abs.(mesh.k) .< mesh.kmax * 2/3 # Dealiasing low-pass filter

        η = zeros(Float64, mesh.N)
        v = zeros(Float64, mesh.N)
		fftη = zeros(Complex{Float64}, mesh.N)
        fftv = zeros(Complex{Float64}, mesh.N)

        new(label, datasize, μ, ϵ, x, F₁, F₂, ∂ₓ, Π⅔, η, v, fftη, fftv )
    end
end


function (m::WhithamBoussinesq)(U::Array{Complex{Float64},2})


    m.fftv .= U[:,2]
	m.fftη .= m.F₂.*U[:,2]
	m.v .= real(ifft(m.fftη))
	m.fftη .= U[:,1]
   	m.η .= real(ifft(U[:,1]))

   	U[:,1] .= -m.∂ₓ.*(m.F₁.*m.fftv.+m.ϵ*m.Π⅔.*m.F₂.*fft(m.η.*m.v))
   	U[:,2] .= -m.∂ₓ.*(m.fftη.+m.ϵ/2*m.Π⅔.*fft(m.v.^2))

end

"""
    mapto(WhithamBoussinesq, data)

"""
function mapto(m::WhithamBoussinesq, data::InitialData)

	[m.Π⅔ .* fft(data.η(m.x)) m.Π⅔ .*fft(data.v(m.x))]

end

"""
    mapfro(WhithamBoussinesq, data)

"""
function mapfro(m::WhithamBoussinesq,
	       datum::Array{Complex{Float64},2})

		   real(ifft(datum[:,1])),real(ifft(datum[:,2]))
end
