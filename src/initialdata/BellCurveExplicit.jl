export BellCurveExplicit

"""
    BellCurveExplicit(param)
    param should contain a value theta

```math
\\eta = 2^{(-|x|^\\theta)}
```

```math
v = 0
```
"""
struct BellCurveExplicit <: InitialData

    η
    v

    function BellCurveExplicit(p :: NamedTuple)
        function η( x :: Vector{Float64} )
            return exp.(.-((abs.(x)).^p.θ).*log(2))
        end
        function v( x :: Vector{Float64} )
            return zeros(length(x))
        end

        new( η , v )

    end

end
