export Random

"""
    Random(param,s,k)

"""
struct Random <: InitialData

    h :: Vector{Float64}
    u :: Vector{Float64}

    function Random(p :: NamedTuple, s :: Int)

        mesh  = Mesh(-p.L, p.L, p.N)


        h=rand(Float64,size(mesh.x))
        u=rand(Float64,size(mesh.x))

        for i in 1:s
            h=cumsum(h.-sum(h)/length(h))
            u=cumsum(u.-sum(u)/length(u))
        end

        h=h/maximum(abs.(h))
        u=u/maximum(abs.(u))

    	new( h,u )

    end


end
