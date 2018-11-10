abstract type InitialData end

export Bump

struct Bump <: InitialData  

    h :: Array{Complex{Float64},1}
    u :: Array{Complex{Float64},1}

    function Bump(p :: Parameters) 

    	mesh  = Mesh(-p.L, p.L, p.N)
    	h = exp.(-(mesh.x).^2)
    	u = zeros(Complex{Float64}, mesh.N)
    	new(h,u)

    end
end
