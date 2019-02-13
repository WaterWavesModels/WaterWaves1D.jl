export CGBSW,mapto,mapfro

"""
    CGBSW( params )

"""
mutable struct CGBSW <: AbstractModel

    label    :: String
    datasize :: Int
    Γ        :: Array{Float64,1}
    Dx       :: Array{Complex{Float64},1}
    H        :: Array{Complex{Float64},1}
    Π⅔       :: BitArray{1}
    ϵ        :: Float64
    hnew     :: Vector{Complex{Float64}}
    unew     :: Vector{Complex{Float64}}

    I₁       :: Vector{Complex{Float64}}
    I₂       :: Vector{Complex{Float64}}
    I₃       :: Vector{Complex{Float64}}

    Px       :: FFTW.FFTWPlan

    function CGBSW( param::NamedTuple)

        label = "Cheng et al."
        datasize = 2
        ϵ = param.ϵ
        mesh  = Mesh(param)
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

        new(label, datasize, Γ, Dx, H, Π⅔, ϵ,
            hnew, unew, I₁, I₂, I₃, Px)

    end
end

function (m::CGBSW)(U::Array{Complex{Float64},2})

    ldiv!(m.hnew, m.Px , view(U,:,1))

    m.I₁  .= view(U,:,2) .* m.Γ
    ldiv!(m.unew, m.Px , m.I₁)
    m.unew .^= 2
    mul!(m.I₁, m.Px , m.unew)
    m.I₁ .*= m.H

    m.I₂  .= view(U,:,1) .* m.Dx
    ldiv!(m.unew, m.Px , m.I₂)
    m.unew .*= m.hnew
    mul!(m.I₂, m.Px , m.unew)

    m.I₃  .= view(U,:,1) .* m.Γ
    ldiv!(m.unew, m.Px, m.I₃)
    m.unew .*= m.hnew
    mul!(m.I₃ , m.Px , m.unew)
    m.I₃ .*= m.H

    m.hnew  .= .- view(U,:,2) .* m.Dx

    m.I₁ .-= m.I₂
    m.I₁ .-= m.I₃
    m.I₁ .*= m.Π⅔
    m.I₁ .*= m.ϵ

    U[:,2] .= view(U,:,1) .* m.H .+ m.I₁
    U[:,1] .= m.hnew

end

"""
    mapto(CGBSW, data)

"""
function mapto(m::CGBSW, data::InitialData)

    [m.Π⅔ .* fft(data.h) m.Π⅔ .* fft(data.u)]

end

"""
    mapfro(CGBSW, data)

"""
function mapfro(m::CGBSW,
           datum ::Array{Complex{Float64},2})

           real(ifft(datum[:,1])),real(ifft(datum[:,2]))
end
