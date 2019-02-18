export BellCurve

"""
    BellCurve(param)
    param should contain a value theta

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

    function BellCurve(p :: NamedTuple)

        mesh  = Mesh(p)
        h     = zeros(Float64, mesh.N)
        h    .= exp.(.-((abs.(mesh.x)).^p.theta).*log(2))
        u     = zeros(Float64, mesh.N)

        new( h, u )

    end

end
