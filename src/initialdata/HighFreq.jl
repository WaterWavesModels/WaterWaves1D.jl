export HighFreq

"""
    HighFreq(param,s,k)

"""
struct HighFreq <: InitialData

    h :: Vector{Float64}
    u :: Vector{Float64}

    function HighFreq(p :: NamedTuple, s :: Real, k :: Int)

    	mesh  = Mesh(-p.L, p.L, p.N)
        h = exp.(-((abs.(mesh.x)).^2)*log(2)).*(cos.(mesh.x).+(1/(2*pi*k)^s)*cos.(2*pi*k*mesh.x))
        u = zeros(Float64, mesh.N)
    	new( h,u )

    end


end
