export HighFreq

"""
    HighFreq(param)
    param should contain s (regularity index) and freq (frequencies)

"""
struct HighFreq <: InitialData

    h :: Vector{Float64}
    u :: Vector{Float64}

    function HighFreq(p :: NamedTuple)

    	mesh  = Mesh(p)
        h = zeros(size(mesh.x))
        for j in p.freq
            if j==0
                h.+= cos.(mesh.x)
            else
                h.+=(1/(2*pi*j)^(p.s))*cos.(2*pi*j*mesh.x)
            end
        end
        h .*= exp.(-((abs.(mesh.x)).^2)*log(2))
        u = mesh.x.*h
    	new( h,u )

    end


end
