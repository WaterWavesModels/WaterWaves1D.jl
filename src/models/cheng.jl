export Cheng

"""
    Cheng( params )

"""
mutable struct Cheng <: AbstractModel

    mesh    :: Mesh
    label   :: String
    data    :: Vector{Tuple{Vector{Complex{Float64}},
			    Vector{Complex{Float64}}}}
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
        
    function Cheng( param::Parameters)

	epsilon = param.ϵ
	mesh  = Mesh(-param.L, param.L, param.N)
        label = "Cheng et al."
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

function (m::Cheng)(h::Vector{Complex{Float64}},
                    u::Vector{Complex{Float64}})
         
    ldiv!(m.hnew, m.Px , h)
    
    m.Int1  .= u
    m.Int1 .*= m.Gamma
    ldiv!(m.unew, m.Px , m.Int1)
    m.unew .^= 2
    mul!(m.Int1, m.Px , m.unew)
    m.Int1 .*= m.H
    
    m.Int2  .= h
    m.Int2 .*= m.Dx
    ldiv!(m.unew, m.Px , m.Int2)
    m.unew .*= m.hnew
    mul!(m.Int2, m.Px , m.unew)
    
    m.Int3  .= h
    m.Int3 .*= m.Gamma
    ldiv!(m.unew, m.Px, m.Int3)
    m.unew .*= m.hnew
    mul!(m.Int3 , m.Px , m.unew)
    m.Int3 .*= m.H
    
    m.hnew  .= -u
    m.hnew .*= m.Dx
    
    m.Int1 .-= m.Int2
    m.Int1 .-= m.Int3
    m.Int1 .*= m.Pi
    m.Int1 .*= m.epsilon
    
    u  .= h
    u .*= m.H
    u .+= m.Int1
    h .= m.hnew
    
end

"""
    init(cheng, data)

"""
function init(m::Cheng, data::InitialData)

    (m.Pi .* fft(data.h), m.Pi .* fft(data.u))

end

"""
    build(cheng, h, u)

"""
function build(m::Cheng, 
	       h::Array{Complex{Float64},1}, 
               u::Array{Complex{Float64},1})

    InitialData(real(ifft(h)),real(ifft(u)))

end
