export WhithamGreenNaghdiKlein,mapto,mapfro

"""
    WhithamGreenNaghdi(params)
	The modified Green-Naghdi model proposed by V. Duchêne, S. Israwi and R. Talhouk.
"""
mutable struct WhithamGreenNaghdiKlein <: AbstractModel

    label   :: String
	datasize:: Int
	μ 		:: Real
	ϵ 		:: Real
	x   	:: Vector{Float64}
    F₀   	:: Vector{Complex{Float64}}
    ∂ₓ      :: Vector{Complex{Float64}}
    Π⅔      :: BitArray{1}
	Id 	    :: BitArray{2}
	FFT 	:: Array{Complex{Float64},2}
	IFFT 	:: Array{Complex{Float64},2}
	IFFTF₀ 	:: Array{Complex{Float64},2}
	M₀      :: Array{Complex{Float64},2}
    h    	:: Vector{Complex{Float64}}
	u    	:: Vector{Complex{Float64}}
	fftv    :: Vector{Complex{Float64}}
	fftη   	:: Vector{Complex{Float64}}
	fftu  	:: Vector{Complex{Float64}}
	I1    	:: Vector{Complex{Float64}}
	I2    	:: Vector{Complex{Float64}}
	L   	:: Array{Complex{Float64},2}
	Precond :: Diagonal{Float64,Array{Float64,1}}
	iterate :: Bool
	ktol 	:: Real
	gtol 	:: Real


    function WhithamGreenNaghdiKlein(param::NamedTuple;iterate=true,SGN=false,dealias=0,ktol=0,gtol=1e-14,precond=true)
		if SGN == true
			label = string("Serre-Green-Naghdi")
		else
			label = string("Whitham-Green-Naghdi")
		end
		datasize = 2
		μ 	= param.μ
		ϵ 	= param.ϵ
		mesh = Mesh(param)
		k = mesh.k
		x 	= mesh.x
		x₀ = mesh.x[1]

		∂ₓ	=  1im * mesh.k
		F₁ 	= tanh.(sqrt(μ)*abs.(k))./(sqrt(μ)*abs.(k))
		F₁[1] 	= 1                 # Differentiation
		if SGN == true
	                F₀ = sqrt(μ)*∂ₓ
	    else
	                F₀ = 1im * sqrt.(3*(1 ./F₁ .- 1)).*sign.(k)
		end
		if precond == true
			Precond = Diagonal( 1 ./  F₁ )
		else
			Precond = Diagonal( 1 .+ μ/3*k.^2 ) #Diagonal( ones(size(k)) )
		end
        Π⅔ 	= abs.(mesh.k) .<= mesh.kmax * (1-dealias/(2+dealias)) # Dealiasing low-pass filter
		FFT = exp.(-1im*k*(x.-x₀)');
        IFFT = exp.(1im*k*(x.-x₀)')/length(x);
		M₀ = IFFT * Diagonal( F₀ ) * FFT
		IFFTF₀ = IFFT * Diagonal( F₀ )
        Id = Diagonal(ones(size(x)));
		h = zeros(Complex{Float64}, mesh.N)
		u, fftv, fftη, fftu, I1, I2 = (similar(h),).*ones(6)
		L = similar(FFT)

        new(label, datasize, μ, ϵ, x, F₀, ∂ₓ, Π⅔, Id, FFT, IFFT, IFFTF₀, M₀, h, u, fftv, fftη, fftu, I1, I2, L, Precond, iterate, ktol, gtol )
    end
end


function (m::WhithamGreenNaghdiKlein)(U::Array{Complex{Float64},2})
	m.fftη .= U[:,1]
	m.h .= 1 .+ m.ϵ*ifft(m.fftη)
	m.fftv .= U[:,2]
	m.fftv[abs.(m.fftv).< m.ktol ].=0   # Krasny filter
	if m.iterate == false
		m.L .= m.Id - 1/3 * m.FFT * Diagonal( 1 ./m.h ) * m.M₀ * Diagonal( m.h.^3 ) * m.IFFTF₀
		m.fftu .= m.L \ m.fftv
	elseif m.iterate == true
        function LL(hatu)
            hatu- 1/3 *fft( 1 ./m.h .* ifft( m.F₀ .* fft( m.h.^3 .* ifft( m.F₀ .* hatu ) ) ) )
		end
		m.fftu .= gmres( LinearMap(LL, length(m.h); issymmetric=false, ismutating=false) , m.fftv ;
				Pl = m.Precond,
				tol = m.gtol )
	end
	m.u .= ifft(m.fftu)
	m.I1 .= m.h.*ifft(m.F₀.*m.fftu)
	m.I2 .= ifft(m.F₀.*fft(m.h.^2 .*m.I1))

   	U[:,1] .= -m.∂ₓ.*fft(m.h .* m.u)
   	U[:,2] .= m.∂ₓ.*(-m.fftη .- m.ϵ/2 * m.Π⅔.*fft(m.u.^2)
					.+ m.ϵ * m.Π⅔.*fft(1/3 * m.u ./ m.h .* m.I2 .+ 1/2 * m.I1.^2))
	#U[:,2] .= -m.∂ₓ.*(m.fftη .+ m.ϵ * m.Π⅔.*fft( m.u.*ifft(m.fftv)
	#				.- 1/2 * m.u.^2 .- 1/2 * m.I1.^2 ) )
	# U[abs.(U).< m.ktol ].=0
end

"""
    mapto(WhithamGreenNaghdiKlein, data)

"""
function mapto(m::WhithamGreenNaghdiKlein, data::InitialData)

	[m.Π⅔ .* fft(data.η(m.x)) m.Π⅔ .*fft(data.v(m.x))]

end

"""
    mapfro(WhithamGreenNaghdi, data)

"""
function mapfro(m::WhithamGreenNaghdiKlein,
	       datum::Array{Complex{Float64},2})
		   m.fftη .= datum[:,1]
	   		m.h .= 1 .+ m.ϵ*ifft(m.fftη)
		   m.L .= m.Id - 1/3 * m.FFT * Diagonal( 1 ./m.h ) * m.M₀ * Diagonal( m.h.^3 ) * m.IFFTF₀

		   real(ifft(datum[:,1])),real(ifft(datum[:,2])),real(ifft(m.L \ datum[:,2]))
end
