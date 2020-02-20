export CGBSW_naive,mapto,mapfro

"""
    CGBSW_naive( params )

"""
mutable struct CGBSW_naive <: AbstractModel

    label   :: String
	datasize:: Int
	x   	:: Array{Float64,1}
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

    function CGBSW_naive( param::NamedTuple)

        label = "Cheng et al. naive"
		datasize = 2
		ϵ = param.ϵ
		mesh  = Mesh(param)
		x = mesh.x
        Γ = abs.(mesh.k)
        Dx    =  1im * mesh.k            # Differentiation
        H     = -1im * sign.(mesh.k)     # Hilbert transform
        Π⅔    = Γ .< mesh.kmax * 2/3 # Dealiasing low-pass filter

        hnew = zeros(Complex{Float64}, mesh.N)
        unew = zeros(Complex{Float64}, mesh.N)

        I₁ = zeros(Complex{Float64}, mesh.N)
        I₂ = zeros(Complex{Float64}, mesh.N)
        I₃ = zeros(Complex{Float64}, mesh.N)

        new(label, datasize, x, Γ, Dx, H, Π⅔, ϵ,
            hnew, unew, I₁, I₂, I₃)

    end
end

function (m::CGBSW_naive)(U::Array{Complex{Float64},2})

	m.hnew .=ifft(U[:,1]);
	m.I₁ .=m.H.*fft(ifft(m.Γ.*U[:,2]).^2);
	m.I₂ .=fft(m.hnew.*ifft(m.Dx.*U[:,1])) ;
	m.I₃ .=m.H.*fft(m.hnew.*ifft(m.Γ.*U[:,1]));
	m.hnew .= U[:,1] ;
	U[:,1] .= -(m.Dx.*U[:,2]) ;
	U[:,2] .= m.H.*m.hnew+m.ϵ*m.Π⅔.*(m.I₁-m.I₂-m.I₃) ;

end

"""
    mapto(CGBSW_naive, data)

"""
function mapto(m::CGBSW_naive, data::InitialData)

	[m.Π⅔ .* fft(data.η(m.x)) m.Π⅔ .*fft(data.v(m.x))]

end

"""
    mapfro(CGBSW_naive, data)

"""
function mapfro(m::CGBSW_naive,
	       datum ::Array{Complex{Float64},2})

		   real(ifft(datum[:,1])),real(ifft(datum[:,2]))
end
