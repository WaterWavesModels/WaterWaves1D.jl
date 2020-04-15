export Matsuno_mod_naive,mapto,mapfro

"""
	Modified Matsuno models with a naive step function
    Matsuno_mod_naive(params)

"""
mutable struct Matsuno_mod_naive <: AbstractModel

    label   :: String
	datasize:: Int
	x   	:: Vector{Float64}
    Γ   	:: Vector{Float64}
    ∂ₓ      :: Vector{Complex{Float64}}
    H       :: Vector{Complex{Float64}}
    Π⅔      :: BitArray{1}
    ϵ 		:: Float64
    hnew    :: Vector{Complex{Float64}}
    unew    :: Vector{Complex{Float64}}
    I₁    :: Vector{Complex{Float64}}
    I₂    :: Vector{Complex{Float64}}
    I₃    :: Vector{Complex{Float64}}

    function Matsuno_mod_naive(param::NamedTuple)

		label = "modified Matsuno"
		datasize = 2
		ϵ 	= param.ϵ
		mesh = Mesh(param)
		x   = mesh.x
        Γ 	= abs.(mesh.k)
    	∂ₓ	=  1im * mesh.k            # Differentiation
        H 	= -1im * sign.(mesh.k)     # Hilbert transform
        Π⅔ 	= Γ .< mesh.kmax * 2/3 # Dealiasing low-pass filter

        hnew = zeros(Complex{Float64}, mesh.N)
        unew = zeros(Complex{Float64}, mesh.N)

        I₁ = zeros(Complex{Float64}, mesh.N)
        I₂ = zeros(Complex{Float64}, mesh.N)
        I₃ = zeros(Complex{Float64}, mesh.N)

        new(label, datasize, x, Γ, ∂ₓ, H, Π⅔, ϵ,
            hnew, unew, I₁, I₂, I₃ )
    end
end


function (m::Matsuno_mod_naive)(U::Array{Complex{Float64},2})

   m.hnew .= ifft(U[:,1])
   m.unew .= ifft(U[:,2])
   m.I₃ .= fft(ifft((m.∂ₓ).*U[:,1]).*exp.(-m.ϵ*ifft((m.Γ).*U[:,1])))
   m.I₁ .= m.H.*U[:,2].-m.ϵ*m.Π⅔.*(m.H.*fft(m.hnew.*ifft(m.Γ.*U[:,2])).+m.∂ₓ.*fft(m.hnew.*m.unew))
   m.I₂ .= -m.Π⅔.*m.I₃-m.ϵ/2*m.Π⅔.*(m.∂ₓ.*fft(m.unew.^2))
   #
   U[:,1] .= m.I₁
   U[:,2] .= m.I₂

end

"""
    mapto(Matsuno, data)
	the velocity should be zero

"""
function mapto(m::Matsuno_mod_naive, data::InitialData)

	[m.Π⅔ .* fft(data.η(m.x)) m.Π⅔ .*fft(data.v(m.x))]

end

"""
    mapfro(Matsuno, data)

"""
function mapfro(m::Matsuno_mod_naive,
	       datum::Array{Complex{Float64},2})

		   real(ifft(datum[:,1])),real(ifft(datum[:,2]))
end
