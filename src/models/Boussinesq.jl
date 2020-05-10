export Boussinesq

"""
    Boussinesq(params)
	params must contain two parameters a and b.
	This computes the abcd Boussinesq model with d=b and c=0. You need a+2*b=1/3 for validity as a long wave model.
"""
mutable struct Boussinesq <: AbstractModel

	label   :: String
	f!		:: Function
	mapto	:: Function
	mapfro	:: Function

    function Boussinesq(param::NamedTuple)

		label = string("Boussinesq")
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

		function f!(U)

		    fftv .= U[:,2]
			fftη .= F₂.*U[:,2]
			v .= real(ifft(fftη))
			fftη .= U[:,1]
		   	η .= real(ifft(U[:,1]))

		   	U[:,1] .= -∂ₓ.*(F₁.*fftv.+ϵ*Π⅔.*F₂.*fft(η.*v))
		   	U[:,2] .= -∂ₓ.*(fftη.+ϵ/2*Π⅔.*fft(v.^2))

		end

		function mapto(data::InitialData)
			[Π⅔ .* fft(data.η(x)) Π⅔ .*fft(data.v(x))]
		end

		function mapfro(U)
			real(ifft(U[:,1])),real(ifft(U[:,2]))
		end

        new(label, f!, mapto, mapfro )
    end
end
