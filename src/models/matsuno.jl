export Matsuno

"""
    Matsuno(params)

"""
mutable struct Matsuno <: AbstractModel
    
    mesh    :: Mesh
    label   :: String
    data    :: Vector{Tuple{Vector{Complex{Float64}},Vector{Complex{Float64}}}}
    Gamma   :: Array{Float64,1}
    Dx      :: Array{Complex{Float64},1}
    H       :: Array{Complex{Float64},1}
    Pi      :: BitArray{1}
    epsilon :: Float64
    hnew    :: Vector{Complex{Float64}}
    unew    :: Vector{Complex{Float64}}
    
    Int1    :: Vector{Complex{Float64}}
    Int2    :: Vector{Complex{Float64}}
    Int3    :: Vector{Complex{Float64}}
    
    Px      :: FFTW.FFTWPlan

    function Matsuno(param::Parameters)

	epsilon = param.ϵ
	mesh  = Mesh(-param.L, param.L, param.N)

        label = "Matsuno"
        data  = []
        Gamma = abs.(mesh.k)
        Dx    =  1im * mesh.k            # Differentiation
        H     = -1im * sign.(mesh.k)     # Hilbert transform
        Pi    = Gamma .< mesh.kmax * 2/3 # Dealiasing low-pass filter

        hnew = zeros(Complex{Float64}, mesh.N)
        unew = zeros(Complex{Float64}, mesh.N)
    
        Int1 = zeros(Complex{Float64}, mesh.N)
        Int2 = zeros(Complex{Float64}, mesh.N)
        Int3 = zeros(Complex{Float64}, mesh.N)
        
        Px  = plan_fft(hnew; flags = FFTW.MEASURE)
      
        new(mesh, label, data, Gamma, Dx, H, Pi, epsilon,
            hnew, unew, Int1, Int2, Int3, Px)
    end
end


function (m::Matsuno)(h::Vector{Complex{Float64}},
                      u::Vector{Complex{Float64}})
    
   # hnew .= real(ifft(h))
   # unew .= real(ifft(u))
   # Int3 .= fft(real(ifft(Dx.*h)).*real(ifft(Gamma.*h)))
   # Int1 .= H.*u-epsilon*Pi.*(H.*fft(hnew.*real(ifft(Gamma.*u))).+Dx.*fft(hnew.*unew))
   # Int2 .= -(Dx.*h)-epsilon/2*Pi.*(Dx.*fft(unew.^2))+epsilon*Pi.*Int3
   # 
   # h .= Int1
   # u .= Int2
     
    m.hnew   .= m.Gamma 
    m.hnew  .*= h         

    ldiv!(m.unew, m.Px, m.hnew )   

    m.hnew   .= m.Dx 
    m.hnew  .*= h            

    ldiv!(m.Int1, m.Px, m.hnew)    

    m.unew  .*= m.Int1            

    mul!(m.Int1, m.Px, m.unew)     

    m.Int1  .*= m.epsilon * m.Pi
    m.Int2   .= m.Dx .* h
    m.Int1  .-= m.Int2

    ldiv!(m.hnew, m.Px, h)
    ldiv!(m.unew, m.Px, u)

    m.Int2   .= m.hnew .* m.unew
    mul!(m.Int3, m.Px, m.Int2)
    m.Int3  .*= m.Dx
    h        .= m.H .* u
    u       .*= m.Gamma
    ldiv!(m.Int2, m.Px, u)
    m.Int2  .*= m.hnew
    mul!(u, m.Px, m.Int2)
    u       .*= m.H
    m.Int3  .+= u
    m.Int3  .*= m.epsilon * m.Pi
    h       .-= m.Int3     
    m.Int3   .= m.unew.^2
    mul!(m.unew, m.Px, m.Int3)
    m.unew  .*= m.Dx
    m.unew  .*= m.epsilon/2 * m.Pi 
    m.Int1  .-= m.unew
    u        .= m.Int1

end

"""
    init( matsuno, data)
"""
function init(m::Matsuno, data::InitialData)

    (m.Pi .* fft(data.h), m.Pi .* fft(data.u))

end

"""
    build( matsuno, h, u)
"""
function build(m::Matsuno, 
	       h::Array{Complex{Float64},1}, 
               u::Array{Complex{Float64},1})

    InitialData(real(ifft(h)),real(ifft(u)))

end
