export Bump

"""
    Bump(param,theta)

```math
h = 2^(-|x|^theta)
```

```math
u = 0
```
"""
struct Bump <: InitialData

    h :: Vector{Float64}
    u :: Vector{Float64}

    function Bump(p :: Parameters,theta :: Real)

    	mesh  = Mesh(-p.L, p.L, p.N)
        h = exp.(-((abs.(mesh.x)).^theta)*log(2))
        u = zeros(Float64, mesh.N)
    	new( h,u )

    end

    function Bump(p :: Parameters)

    	mesh  = Mesh(-p.L, p.L, p.N)
        h = exp.(-(mesh.x).^2)
        u = zeros(Float64, mesh.N)
    	new( h,u )

    end

end
