export BellCurve

"""
    BellCurve(param,theta)

```math
h = 2^{(-|x|^\\theta)}
```

```math
u = 0
```
"""
struct BellCurve <: InitialData

    h :: Vector{Float64}
    u :: Vector{Float64}

    function BellCurve(p :: NamedTuple,theta :: Real)

    	mesh  = Mesh(-p.L, p.L, p.N)
        h = zeros(Float64, mesh.N)
        h .= exp.(.-((abs.(mesh.x)).^theta).*log(2))
        u = zeros(Float64, mesh.N)
    	new( h,u )

    end

    function BellCurve(p :: NamedTuple)

    	mesh  = Mesh(-p.L, p.L, p.N)
        h = zeros(Float64, mesh.N)
        h .= exp.(.-(mesh.x).^2)
        u = zeros(Float64, mesh.N)
    	new( h,u )

    end

end
