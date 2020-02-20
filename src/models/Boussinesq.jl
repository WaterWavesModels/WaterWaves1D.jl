export Boussinesq,mapto,mapfro

"""
    Boussinesq(params)
	A Hamiltonian Boussinesq system derived in BonaSmith76 and studied in BonaChenSaut02

"""
mutable struct Boussinesq <: AbstractModel

    label   :: String
	datasize:: Int
	μ 		:: Float64
	ϵ 		:: Float64
	x   	:: Array{Float64,1}
    G₁   	:: Array{Float64,1}
    ∂ₓ      :: Array{Complex{Float64},1}
    Π⅔      :: BitArray{1}
	η    	:: Vector{Float64}
	v    	:: Vector{Float64}
    fftη    :: Vector{Complex{Float64}}
    fftv    :: Vector{Complex{Float64}}

    function Boussinesq(param::NamedTuple)

		label = "Boussinesq"
		datasize = 2
		μ 	= param.μ
		ϵ 	= param.ϵ
		mesh = Mesh(param)
		x = mesh.x
        G₁ 	= 1 ./(1 .+μ/3*abs.(mesh.k).^2)
    	∂ₓ	=  1im * mesh.k            # Differentiation
        Π⅔ 	= abs.(mesh.k) .< mesh.kmax * 2/3 # Dealiasing low-pass filter

        η = zeros(Float64, mesh.N)
        v = zeros(Float64, mesh.N)
		fftη = zeros(Complex{Float64}, mesh.N)
        fftv = zeros(Complex{Float64}, mesh.N)

        new(label, datasize, μ, ϵ, x, G₁, ∂ₓ, Π⅔, η, v, fftη, fftv )
    end
end


function (m::Boussinesq)(U::Array{Complex{Float64},2})


	m.fftη .= U[:,1]
    m.fftv .= U[:,2]
   	m.η .= real(ifft(U[:,1]))
   	m.v .= real(ifft(U[:,2]))

   	U[:,1] .= -m.∂ₓ.*(m.fftv.+m.ϵ*m.Π⅔.*m.G₁.*fft(m.η.*m.v))
   	U[:,2] .= -m.∂ₓ.*m.G₁.*(m.fftη.+m.ϵ/2*m.Π⅔.*fft(m.v.^2))

end

"""
    mapto(Boussinesq, data)

"""
function mapto(m::Boussinesq, data::InitialData)

	[m.Π⅔ .* fft(data.η(m.x)) m.Π⅔ .*m.G₁.*fft(data.v(m.x))]

end

"""
    mapfro(Boussinesq, data)

"""
function mapfro(m::Boussinesq,
	       datum::Array{Complex{Float64},2})
		   G = m.G₁.^(-1)
		   real(ifft(datum[:,1])),real(ifft(G.* datum[:,2]))
end
