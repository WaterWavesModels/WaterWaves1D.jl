export WhithamBoussinesq

"""
    WhithamBoussinesq(params)
	params must contain a parameter α.
	If α = 1, then the model has been introduced and studied by E. Dinvay and collaborators.
	If α = 1/2, then the model is a quasilinear version.
	If α < 1/2, then expect instabilities stemming from ill-posedness of the model.
"""
mutable struct WhithamBoussinesq <: AbstractModel

	label   :: String
	f!		:: Function
	mapto	:: Function
	mapfro	:: Function

    function WhithamBoussinesq(param::NamedTuple)

		label = string("Whitham-Boussinesq")
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


        new(label, f!, mapto, mapfro)
    end
end
