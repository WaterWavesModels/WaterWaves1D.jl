export WhithamGreenNaghdiSym,mapto,mapfro

"""
    WhithamGreenNaghdiSym(params)
	The modified Green-Naghdi model proposed by V. Duchêne, S. Israwi and R. Talhouk.
"""
mutable struct WhithamGreenNaghdiSym <: AbstractModel

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
	F₀FFT 	:: Array{Complex{Float64},2}
    h    	:: Vector{Complex{Float64}}
	u    	:: Vector{Complex{Float64}}
	v    	:: Vector{Complex{Float64}}
	fftη    :: Vector{Complex{Float64}}
	fftu   :: Vector{Complex{Float64}}
	ffthv  	:: Vector{Complex{Float64}}
	hdu    	:: Vector{Complex{Float64}}
	L   	:: Array{Complex{Float64},2}
	Precond :: Diagonal{Float64,Array{Float64,1}}
	iterate :: Bool
	ktol 	:: Real
	gtol 	:: Real


    function WhithamGreenNaghdiSym(param::NamedTuple;iterate=true,SGN=false,dealias=0,ktol=0,gtol=1e-14,precond=true)
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
			Precond = Diagonal( 1 .+ μ*k.^2 ) #Diagonal( ones(size(k)) )
		end
        Π⅔ 	= abs.(mesh.k) .<= mesh.kmax * (1-dealias/(2+dealias)) # Dealiasing low-pass filter
		FFT = exp.(-1im*k*(x.-x₀)')
        IFFT = exp.(1im*k*(x.-x₀)')/length(x)
		IFFTF₀ = IFFT * Diagonal( F₀ )
		F₀FFT = Diagonal( F₀ ) * FFT
        Id = Diagonal(ones(size(x)));
		h = zeros(Complex{Float64}, mesh.N)
		u, v, fftη, fftu, ffthv, hdu = (similar(h),).*ones(6)
		L = similar(FFT)

        new(label, datasize, μ, ϵ, x, F₀, ∂ₓ, Π⅔, Id, FFT, IFFT, IFFTF₀, F₀FFT, h, u, v, fftη, fftu, ffthv, hdu, L, Precond, iterate, ktol, gtol )
    end
end


function (m::WhithamGreenNaghdiSym)(U::Array{Complex{Float64},2})
	m.fftη .= U[:,1]
	m.v .= ifft(U[:,2])
	#m.fftη[abs.(m.fftη).< m.ktol ].=0   # Krasny filter
	#m.fftv[abs.(m.fftv).< m.ktol ].=0   # Krasny filter
	m.h .= 1 .+ m.ϵ*ifft(m.fftη)
	m.ffthv .= fft(m.h.*m.v)
	if m.iterate == false
		m.L .= m.FFT * Diagonal( m.h ) * m.IFFT - 1/3 * m.F₀FFT * Diagonal( m.h.^3 ) * m.IFFTF₀
		m.fftu .= m.L \ m.ffthv
	elseif m.iterate == true
        function LL(hatu)
            fft( m.h .* ifft(hatu) )- 1/3 * m.F₀ .* fft( m.h.^3 .* ifft( m.F₀ .* hatu ) )
		end
		m.fftu .= gmres( LinearMap(LL, length(m.h); issymmetric=false, ismutating=false) , m.ffthv ;
				Pl = m.Precond,
				tol = m.gtol )
	end
	m.u .= ifft(m.fftu)
	m.hdu .= m.h .* ifft(m.F₀.*m.fftu)

   	U[:,1] .= -m.∂ₓ.*fft(m.h.*m.u)
   	U[:,2] .= -m.∂ₓ.*(m.fftη .+ m.ϵ * m.Π⅔.*fft( m.u.*m.v
					.- 1/2 * m.u.^2 .- 1/2 * m.hdu.^2 ) )
	# U[abs.(U).< m.ktol ].=0 # Krasny filter
end

"""
    mapto(WhithamGreenNaghdiSym, data)

"""
function mapto(m::WhithamGreenNaghdiSym, data::InitialData)

	[m.Π⅔ .* fft(data.η(m.x)) m.Π⅔ .*fft(data.v(m.x))]

end

"""
    mapfro(WhithamGreenNaghdiSym, data)

"""
function mapfro(m::WhithamGreenNaghdiSym,
	       datum::Array{Complex{Float64},2})
	   		m.h .= 1 .+ m.ϵ*ifft(datum[:,1])
		    m.L .= m.Id - 1/3 * m.FFT * Diagonal( 1 ./m.h ) * m.IFFT*  m.F₀FFT * Diagonal( m.h.^3 ) * m.IFFTF₀

		   real(ifft(datum[:,1])),real(ifft(datum[:,2])),real(ifft(m.L \ datum[:,2]))
end
