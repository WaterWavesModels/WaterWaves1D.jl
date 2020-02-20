export Init

"""
    Init(init)
    init is a NamedTuple and should contain a value η and a value v

```math
\\eta = 2^{(-|x|^\\theta)}
```

```math
v = 0
```
"""
struct Init <: InitialData

    η
    v

    function Init(p :: NamedTuple)
        new( x->p.η(x) , x->p.v(x) )
    end

    function Init(η , v)
        new( x->η(x) , x->v(x) )
    end

end
