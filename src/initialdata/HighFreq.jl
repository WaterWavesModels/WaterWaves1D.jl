export HighFreq

"""
    HighFreq(param)
    param should contain s (regularity index) and freq (frequencies)

"""
struct HighFreq <: InitialData

    η
    v

    function HighFreq(p :: NamedTuple)

        function η( x :: Vector{Float64} )

            h = zeros(size(x))
            for j in p.freq
                if j==0
                    h.+= cos.(x)
                else
                    h.+=(1/(2*pi*j)^(p.s))*cos.(2*pi*j*x)
                end
            end
            h .*= exp.(-((abs.(x)).^2)*log(2))
            return h
        end

        function v( x :: Vector{Float64} )
            return u = x.*η(x)
        end
    	new( η,v )

    end


end
