using FFTW, LinearAlgebra

export CGBSW,construct,reconstruct

"""
    CGBSW( params )

"""
mutable struct CGBSW <: AbstractModel

    mesh    :: Mesh
    label   :: String
    data    :: Vector{Tuple{Vector{Complex{Float64}},
			    Vector{Complex{Float64}}}}
    Γ   	:: Array{Float64,1}
    Dx      :: Array{Complex{Float64},1}
    H       :: Array{Complex{Float64},1}
    Π⅔      :: BitArray{1}
    ϵ 		:: Float64
  	hnew 	:: Vector{Complex{Float64}}
    unew    :: Vector{Complex{Float64}}

    I₁   	:: Vector{Complex{Float64}}
    I₂    	:: Vector{Complex{Float64}}
    I₃    	:: Vector{Complex{Float64}}

    Px      :: FFTW.FFTWPlan

    function CGBSW( param::NamedTuple)

	ϵ = param.ϵ
	mesh  = Mesh(-param.L, param.L, param.N)
        label = "Cheng et al."
        data  = []
        Γ = abs.(mesh.k)
        Dx    =  1im * mesh.k            # Differentiation
        H     = -1im * sign.(mesh.k)     # Hilbert transform
        Π⅔    = Γ .< mesh.kmax * 2/3 # Dealiasing low-pass filter

        hnew = zeros(Complex{Float64}, mesh.N)
        unew = zeros(Complex{Float64}, mesh.N)

        I₁ = zeros(Complex{Float64}, mesh.N)
        I₂ = zeros(Complex{Float64}, mesh.N)
        I₃ = zeros(Complex{Float64}, mesh.N)

        Px  = plan_fft(hnew; flags = FFTW.MEASURE)

        new(mesh, label, data, Γ, Dx, H, Π⅔, ϵ,
            hnew, unew, I₁, I₂, I₃, Px)

    end
end

function (m::CGBSW)(h::Vector{Complex{Float64}},
                    u::Vector{Complex{Float64}})

    ldiv!(m.hnew, m.Px , h)

    m.I₁  .= u
    m.I₁ .*= m.Γ
    ldiv!(m.unew, m.Px , m.I₁)
    m.unew .^= 2
    mul!(m.I₁, m.Px , m.unew)
    m.I₁ .*= m.H

    m.I₂  .= h
    m.I₂ .*= m.Dx
    ldiv!(m.unew, m.Px , m.I₂)
    m.unew .*= m.hnew
    mul!(m.I₂, m.Px , m.unew)

    m.I₃  .= h
    m.I₃ .*= m.Γ
    ldiv!(m.unew, m.Px, m.I₃)
    m.unew .*= m.hnew
    mul!(m.I₃ , m.Px , m.unew)
    m.I₃ .*= m.H

    m.hnew  .= -u
    m.hnew .*= m.Dx

    m.I₁ .-= m.I₂
    m.I₁ .-= m.I₃
    m.I₁ .*= m.Π⅔
    m.I₁ .*= m.ϵ

    u  .= h
    u .*= m.H
    u .+= m.I₁
    h .= m.hnew

end

"""
    construct(CGBSW, data)

"""
function construct(m::CGBSW, data::InitialData)

    (m.Π⅔ .* fft(data.h), m.Π⅔ .* fft(data.u))

end

"""
    reconstruct(CGBSW, h, u)

"""
function reconstruct(m::CGBSW,
	       h::Array{Complex{Float64},1},
               u::Array{Complex{Float64},1})

    InitialData(real(ifft(h)),real(ifft(u)))

end
