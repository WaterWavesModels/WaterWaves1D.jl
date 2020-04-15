export Random

"""
    Random(param)

param should contain s (regularity index)

"""
struct Random <: InitialData

    η
    v

    function Random(p :: NamedTuple)

        function η( x :: Vector{Float64} )
            h = rand(Float64,size(x))

            for i in 1:p.s
                h .= cumsum(h.-sum(h)/length(h))
            end

            return h ./ maximum(abs.(h))
        end
        function v( x :: Vector{Float64} )
            h = rand(Float64,size(x))

            for i in 1:p.s
                h .= cumsum(h.-sum(h)/length(h))
            end

            return h ./ maximum(abs.(h))
        end

    	new( η,v )

    end

end
