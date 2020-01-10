export fdBoussinesq_2,mapto,mapfro

"""
    fdBoussinesq_1(params)
	A quasilinear version of the model introduced and studied by E. Dinvay and collaborators

"""
mutable struct fdBoussinesq_2 <: AbstractModel

    label   :: String
	datasize:: Int
    G₁   	:: Array{Float64,1}
	G₂   	:: Array{Float64,1}
    ∂ₓ      :: Array{Complex{Float64},1}
    Π⅔      :: BitArray{1}
	μ 		:: Float64
    ϵ 		:: Float64
	h    	:: Vector{Float64}
	v    	:: Vector{Float64}
    ffth    :: Vector{Complex{Float64}}
    fftv    :: Vector{Complex{Float64}}

    function fdBoussinesq_2(param::NamedTuple)

		label = "fdBoussinesq_2"
		datasize = 2
		μ 	= param.μ
		ϵ 	= param.ϵ
		mesh = Mesh(param)
		#G₀ 	= tanh(sqrt(μ)*abs.(mesh.k)).*(sqrt(μ)*abs.(mesh.k))
        G₁ 	= sqrt.(tanh.(sqrt(μ)*abs.(mesh.k))./(sqrt(μ)*abs.(mesh.k)))
		G₁[1] 	= 1
		G₂ 	= G₁.^2
    	∂ₓ	=  1im * mesh.k            # Differentiation
        Π⅔ 	= abs.(mesh.k) .< mesh.kmax * 2/3 # Dealiasing low-pass filter

		h = zeros(Float64, mesh.N)
        v = zeros(Float64, mesh.N)
		ffth = zeros(Complex{Float64}, mesh.N)
        fftv = zeros(Complex{Float64}, mesh.N)

        new(label, datasize, G₁, G₂, ∂ₓ, Π⅔, μ, ϵ, h, v, ffth, fftv )
    end
end


function (m::fdBoussinesq_2)(U::Array{Complex{Float64},2})

	m.ffth .= U[:,1]
    m.fftv .= U[:,2]
   	m.h .= real(ifft(U[:,1]))
   	m.v .= real(ifft(m.G₁.*U[:,2]))

   U[:,1] .= -m.∂ₓ.*(m.G₂.*m.fftv.+m.ϵ*m.Π⅔.*m.G₁.*fft(m.h.*m.v))
   U[:,2] .= -m.∂ₓ.*(m.ffth.+m.ϵ/2*m.Π⅔.*fft(m.v.^2))

end

"""
    mapto(fdBoussinesq_2, data)

"""
function mapto(m::fdBoussinesq_2, data::InitialData)

    [m.Π⅔ .* fft(data.h) m.Π⅔ .* fft(data.u)]

end

"""
    mapfro(fdBoussinesq_2, data)

"""
function mapfro(m::fdBoussinesq_2,
	       datum::Array{Complex{Float64},2})

		   real(ifft(datum[:,1])),real(ifft(datum[:,2]))
end
