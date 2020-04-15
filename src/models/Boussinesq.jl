export Boussinesq,mapto,mapfro

"""
    Boussinesq(params)
	params must contain two parameters a and b.
	This computes the abcd Boussinesq model with d=b and c=0. You need a+2*b=1/3 for validity as a long wave model.
"""
mutable struct Boussinesq <: AbstractModel

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

    function Boussinesq(param::NamedTuple)

		label = string("Boussinesq")
		datasize = 2
		μ 	= param.μ
		ϵ 	= param.ϵ
		mesh = Mesh(param)
		x 	= mesh.x
		F₂ = 1 ./(1 .+μ*param.b*abs.(mesh.k).^2)
		F₁ 	= (1 .-μ*param.a*abs.(mesh.k).^2).*(F₂.^2)
    	∂ₓ	=  1im * mesh.k            # Differentiation
        Π⅔ 	= abs.(mesh.k) .< mesh.kmax * 2/3 # Dealiasing low-pass filter

        η = zeros(Float64, mesh.N)
        v = zeros(Float64, mesh.N)
		fftη = zeros(Complex{Float64}, mesh.N)
        fftv = zeros(Complex{Float64}, mesh.N)

        new(label, datasize, μ, ϵ, x, F₁, F₂, ∂ₓ, Π⅔, η, v, fftη, fftv )
    end
end


function (m::Boussinesq)(U::Array{Complex{Float64},2})


    m.fftv .= U[:,2]
	m.fftη .= m.F₂.*U[:,2]
	m.v .= real(ifft(m.fftη))
	m.fftη .= U[:,1]
   	m.η .= real(ifft(U[:,1]))

   	U[:,1] .= -m.∂ₓ.*(m.F₁.*m.fftv.+m.ϵ*m.Π⅔.*m.F₂.*fft(m.η.*m.v))
   	U[:,2] .= -m.∂ₓ.*(m.fftη.+m.ϵ/2*m.Π⅔.*fft(m.v.^2))

end

"""
    mapto(Boussinesq, data)

"""
function mapto(m::Boussinesq, data::InitialData)

	[m.Π⅔ .* fft(data.η(m.x)) m.Π⅔ .*fft(data.v(m.x))]

end

"""
    mapfro(Boussinesq, data)

"""
function mapfro(m::Boussinesq,
	       datum::Array{Complex{Float64},2})

		   real(ifft(datum[:,1])),real(ifft(datum[:,2]))
end
