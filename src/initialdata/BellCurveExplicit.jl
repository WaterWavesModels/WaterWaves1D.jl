export BellCurveExplicit

"""
    BellCurveExplicit(param)
    param should contain a value theta

```math
h = 2^{(-|x|^\\theta)}
```

```math
u = 0
```
"""
struct BellCurveExplicit <: InitialData

    init

    function BellCurveExplicit(p :: NamedTuple)
        θ = p.theta
        function init( x :: Vector{Float64} )
            return exp.(.-((abs.(x)).^θ).*log(2))
        end

        new( init )

    end

end
