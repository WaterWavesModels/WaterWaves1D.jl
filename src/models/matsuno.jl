export Matsuno,mapto,mapfro

"""
    Matsuno(params)

"""
mutable struct Matsuno <: AbstractModel

    mesh    :: Mesh
    label   :: String
    Γ   	:: Array{Float64,1}
    Dx      :: Array{Complex{Float64},1}
    H       :: Array{Complex{Float64},1}
    Π⅔      :: BitArray{1}
    ϵ 		:: Float64
    hnew    :: Vector{Complex{Float64}}
    unew    :: Vector{Complex{Float64}}
	I₀    :: Vector{Complex{Float64}}
    I₁    :: Vector{Complex{Float64}}
    I₂    :: Vector{Complex{Float64}}
    I₃    :: Vector{Complex{Float64}}

    Px      :: FFTW.FFTWPlan

    function Matsuno(param::NamedTuple)

		ϵ = param.ϵ
		mesh  = Mesh(-param.L, param.L, param.N)
        label = "Matsuno"
        Γ = abs.(mesh.k)
        Dx    =  1im * mesh.k            # Differentiation
        H     = -1im * sign.(mesh.k)     # Hilbert transform
        Π⅔    = Γ .< mesh.kmax * 2/3 # Dealiasing low-pass filter

        hnew = zeros(Complex{Float64}, mesh.N)
        unew = zeros(Complex{Float64}, mesh.N)

		I₀ = zeros(Complex{Float64}, mesh.N)
        I₁ = zeros(Complex{Float64}, mesh.N)
        I₂ = zeros(Complex{Float64}, mesh.N)
        I₃ = zeros(Complex{Float64}, mesh.N)

        Px  = plan_fft(hnew; flags = FFTW.MEASURE)

        new(mesh, label, Γ, Dx, H, Π⅔, ϵ,
            hnew, unew, I₀, I₁, I₂, I₃, Px)
    end
end


function (m::Matsuno)(U::Array{Complex{Float64},2})

   # hnew .= real(ifft(h))
   # unew .= real(ifft(u))
   # I₃ .= fft(real(ifft(Dx.*h)).*real(ifft(Γ.*h)))
   # I₁ .= H.*u-ϵ*Π⅔.*(H.*fft(hnew.*real(ifft(Γ.*u))).+Dx.*fft(hnew.*unew))
   # I₂ .= -(Dx.*h)-ϵ/2*Π⅔.*(Dx.*fft(unew.^2))+ϵ*Π⅔.*I₃
   #
   # h .= I₁
   # u .= I₂

    m.hnew   .= m.Γ
    m.hnew  .*= U[:,1]

    ldiv!(m.unew, m.Px, m.hnew )

    m.hnew   .= m.Dx
    m.hnew  .*= U[:,1]

    ldiv!(m.I₁, m.Px, m.hnew)

    m.unew  .*= m.I₁

    mul!(m.I₁, m.Px, m.unew)

    m.I₁  .*= m.ϵ * m.Π⅔
    m.I₂   .= m.Dx .* U[:,1]
    m.I₁  .-= m.I₂

    ldiv!(m.hnew, m.Px, U[:,1])
    ldiv!(m.unew, m.Px, U[:,2])

    m.I₂   .= m.hnew .* m.unew
    mul!(m.I₃, m.Px, m.I₂)
    m.I₃  .*= m.Dx
    U[:,1]  .= m.H .* U[:,2]
	m.I₀ .= m.Γ.*U[:,2]
    ldiv!(m.I₂, m.Px, m.I₀)
    m.I₂  .*= m.hnew
    mul!(m.hnew, m.Px, m.I₂)
    m.hnew  .*= m.H
    m.I₃  .+= m.hnew
    m.I₃  .*= m.ϵ * m.Π⅔
    U[:,1] .-= m.I₃
    m.I₃   .= m.unew.^2
    mul!(m.unew, m.Px, m.I₃)
    m.unew  .*= m.Dx
    m.unew  .*= m.ϵ/2 * m.Π⅔
    m.I₁  .-= m.unew
    U[:,2]  .= m.I₁

end

"""
    mapto(Matsuno, data)

"""
function mapto(m::Matsuno, data::InitialData)

    [m.Π⅔ .* fft(data.h) m.Π⅔ .* fft(data.u)]

end

"""
    mapfro(Matsuno, data)

"""
function mapfro(m::Matsuno,
	       data::Array{Complex{Float64},2})

		   real(ifft(data[:,1])),real(ifft(data[:,2]))
end
