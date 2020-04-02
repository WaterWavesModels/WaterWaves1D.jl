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
	M₀  	:: Array{Complex{Float64},2}
    h    	:: Vector{Complex{Float64}}
	u    	:: Vector{Complex{Float64}}
	v    	:: Vector{Complex{Float64}}
	η       :: Vector{Complex{Float64}}
	fftu    :: Vector{Complex{Float64}}
	ffthv  	:: Vector{Complex{Float64}}
	hdu    	:: Vector{Complex{Float64}}
	L   	:: Array{Complex{Float64},2}
	Precond :: Any
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
		if SGN == true
			F₁ 	= 1 ./ (1 .+ μ/3*k.^2)
			F₀ = sqrt(μ)*∂ₓ
	    else
			F₁ 	= tanh.(sqrt(μ)*abs.(k))./(sqrt(μ)*abs.(k))
			F₁[1] 	= 1
			F₀ = 1im * sqrt.(3*(1 ./F₁ .- 1)).*sign.(k)
		end
        Π⅔ 	= abs.(mesh.k) .<= mesh.kmax * (1-dealias/(2+dealias)) # Dealiasing low-pass filter
		FFT = exp.(-1im*k*(x.-x₀)')
        IFFT = exp.(1im*k*(x.-x₀)')/length(x)
		M₀ = IFFT * Diagonal( F₀ ) *IFFT
        Id = Diagonal(ones(size(x)));
		if precond == true
			Precond = lu(IFFT*Diagonal( 1 ./ F₁ )*FFT)
		else
			Precond = Identity() #Diagonal( ones(size(k)) )
		end
		h = zeros(Complex{Float64}, mesh.N)
		u, v, η, fftu, ffthv, hdu = (copy(h),).*ones(6)
		L = similar(FFT)

        new(label, datasize, μ, ϵ, x, F₀, ∂ₓ, Π⅔, Id, M₀, h, u, v, η, fftu, ffthv, hdu, L, Precond, iterate, ktol, gtol )
    end
end


function (m::WhithamGreenNaghdiSym)(U::Array{Complex{Float64},2})
	m.η .= U[:,1]
	m.v .= U[:,2]
	#m.η[abs.(m.η).< m.ktol ].=0   # Krasny filter
	#m.fftv[abs.(m.fftv).< m.ktol ].=0   # Krasny filter
	m.h .= 1 .+ m.ϵ*U[:,1]
	#m.ffthv .= fft(m.h.*m.v)
	if m.iterate == false
		m.L .= Symmetric(Diagonal( m.h ) - 1/3 * m.M₀ * Diagonal( m.h.^3 ) * m.M₀)
		m.u .= m.L \ (m.h.*m.v)
	elseif m.iterate == true
        function LL(u)
            m.h .* u - 1/3 * ifft( m.F₀ .* fft( m.h.^3 .* ifft( m.F₀ .* fft( u ) ) ) )
		end
		#m.u = m.v
		cg!( m.u, LinearMap(LL, length(m.h); issymmetric=true) , (m.h.*m.v) ;
				Pl = m.Precond,
				verbose = false,
				tol = m.gtol )
	end
	#m.u .= ifft(m.fftu)
	m.hdu .= m.h .* ifft(m.F₀.* fft(m.u))

   	U[:,1] .= -ifft(m.∂ₓ.*m.Π⅔.*fft(m.h.*m.u))
   	U[:,2] .= -ifft(m.∂ₓ.*m.Π⅔.*fft(m.η .+ m.ϵ * m.u.*m.v
					.- m.ϵ/2 * m.u.^2 .- m.ϵ/2 * m.hdu.^2 ) )
	# U[abs.(U).< m.ktol ].=0 # Krasny filter
end

"""
    mapto(WhithamGreenNaghdiSym, data)

"""
function mapto(m::WhithamGreenNaghdiSym, data::InitialData)
	m.η .= Complex.(data.η(m.x))
	m.v .= Complex.(data.v(m.x))
	m.h .= 1 .+ m.ϵ*m.η
	m.L .= Symmetric(Diagonal( m.h ) - 1/3 * m.M₀ * Diagonal( m.h.^3 ) * m.M₀)
	m.u .= m.L \ (m.h.*m.v)
	return [m.η m.v]

end

"""
    mapfro(WhithamGreenNaghdiSym, data)

"""
function mapfro(m::WhithamGreenNaghdiSym,
	       datum::Array{Complex{Float64},2})
	   		m.h .= 1 .+ m.ϵ*datum[:,1]
			m.L .= Symmetric(Diagonal( m.h ) - 1/3 * m.M₀ * Diagonal( m.h.^3 ) * m.M₀)

		   real.(datum[:,1]),real.(datum[:,2]),real.(m.L \ (m.h.*datum[:,2]))
end
