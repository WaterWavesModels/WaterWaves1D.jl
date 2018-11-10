abstract type InitialData end

export Bump

"""
    Bump(param)

```math
h = exp(-x^2)
```

```math
u = 0
```
"""
struct Bump <: InitialData  

    mesh :: Mesh

    function Bump(p :: Parameters) 

    	mesh  = Mesh(-p.L, p.L, p.N)
    	new( mesh )

    end
end


function(b::Bump)(h, u)

    h .= exp.(-(b.mesh.x).^2)
    u .= zeros(Complex{Float64}, b.mesh.N)

end

