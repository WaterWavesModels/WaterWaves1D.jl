using FFTW, LinearAlgebra

export Matsuno,construct,reconstruct

"""
    Matsuno(params)

"""
mutable struct Matsuno <: AbstractModel

    mesh    :: Mesh
    label   :: String
    data    :: Vector{Tuple{Vector{Complex{Float64}},Vector{Complex{Float64}}}}
    Γ   	:: Array{Float64,1}
    Dx      :: Array{Complex{Float64},1}
    H       :: Array{Complex{Float64},1}
    Π⅔      :: BitArray{1}
    ϵ 		:: Float64
    hnew    :: Vector{Complex{Float64}}
    unew    :: Vector{Complex{Float64}}

    I₁    :: Vector{Complex{Float64}}
    I₂    :: Vector{Complex{Float64}}
    I₃    :: Vector{Complex{Float64}}

    Px      :: FFTW.FFTWPlan

    function Matsuno(param::NamedTuple)

	ϵ = param.ϵ
	mesh  = Mesh(-param.L, param.L, param.N)

        label = "Matsuno"
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


function (m::Matsuno)(h::Vector{Complex{Float64}},
                      u::Vector{Complex{Float64}})

   # hnew .= real(ifft(h))
   # unew .= real(ifft(u))
   # I₃ .= fft(real(ifft(Dx.*h)).*real(ifft(Γ.*h)))
   # I₁ .= H.*u-ϵ*Π⅔.*(H.*fft(hnew.*real(ifft(Γ.*u))).+Dx.*fft(hnew.*unew))
   # I₂ .= -(Dx.*h)-ϵ/2*Π⅔.*(Dx.*fft(unew.^2))+ϵ*Π⅔.*I₃
   #
   # h .= I₁
   # u .= I₂

    m.hnew   .= m.Γ
    m.hnew  .*= h

    ldiv!(m.unew, m.Px, m.hnew )

    m.hnew   .= m.Dx
    m.hnew  .*= h

    ldiv!(m.I₁, m.Px, m.hnew)

    m.unew  .*= m.I₁

    mul!(m.I₁, m.Px, m.unew)

    m.I₁  .*= m.ϵ * m.Π⅔
    m.I₂   .= m.Dx .* h
    m.I₁  .-= m.I₂

    ldiv!(m.hnew, m.Px, h)
    ldiv!(m.unew, m.Px, u)

    m.I₂   .= m.hnew .* m.unew
    mul!(m.I₃, m.Px, m.I₂)
    m.I₃  .*= m.Dx
    h        .= m.H .* u
    u       .*= m.Γ
    ldiv!(m.I₂, m.Px, u)
    m.I₂  .*= m.hnew
    mul!(u, m.Px, m.I₂)
    u       .*= m.H
    m.I₃  .+= u
    m.I₃  .*= m.ϵ * m.Π⅔
    h       .-= m.I₃
    m.I₃   .= m.unew.^2
    mul!(m.unew, m.Px, m.I₃)
    m.unew  .*= m.Dx
    m.unew  .*= m.ϵ/2 * m.Π⅔
    m.I₁  .-= m.unew
    u        .= m.I₁

end

"""
    construct( matsuno, data)
"""
function construct(m::Matsuno, data::InitialData)

    (m.Π⅔ .* fft(data.h), m.Π⅔ .* fft(data.u))

end

"""
    reconstruct( matsuno, h, u)
"""
function reconstruct(m::Matsuno,
	       h::Array{Complex{Float64},1},
               u::Array{Complex{Float64},1})

    InitialData(real(ifft(h)),real(ifft(u)))

end
