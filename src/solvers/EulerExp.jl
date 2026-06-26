export EulerExp, EulerExp_naive
export step!

"""
    EulerExp(arguments;realdata)

Explicit Exponential Euler solver.

Construct an object of type `TimeSolver` to be used in `Problem(model, initial, param; solver::TimeSolver)`

Arguments can be either
0. an object of type `AbstractModel`, which contains the necessary;
1. an `Array` of size `(N,datasize)` where `N` is the number of collocation points and `datasize` the number of equations solved;
2. `(param,datasize)` where `param is a `NamedTuple` containing a key `N`, and `datasize` a integer (optional, by default `datasize=2`).

The keyword argument `realdata` is optional, and determines whether pre-allocated vectors are real- or complex-valued.
By default, they are either determined by the model or the type of the array in case `0.` and `1.`, complex-valued in case `2.`.

"""
struct EulerExp <: TimeSolver

    U1    :: Array
    φ     :: Function
    label :: String

    function EulerExp( U :: Array; realdata=nothing )
        U1 = deepcopy(U)
        φ(z)=(exp(z+eps())-1)/(z+eps())
        if realdata==true
            U1 = real.(U1);
        end
        if realdata==false
            U1 = complex.(U1);
        end
        new( U1, φ, "exponential Euler" )
    end

    function EulerExp( model :: AbstractModel; realdata=nothing )
        U=model.mapto(Init(x->0*x,x->0*x))
        EulerExp( U; realdata=realdata)
    end
    function EulerExp( param::NamedTuple, systemsize=2::Int; realdata=nothing )
        EulerExp( [Array{Complex{Float64}}(undef,param.N) for _ in 1:systemsize]  ; realdata=realdata)
    end
    function EulerExp( datasize, systemsize=2::Int; realdata=nothing )
        EulerExp( [Array{Complex{Float64}}(undef,datasize) for _ in 1:systemsize] ; realdata=realdata)
    end
end

@doc raw"""
    step!(solver :: EulerExp, model :: AbstractModel , U, dt)

Integration step of the exponential Euler solver applied to solutions to the equation `u'=A u + g(u)`.

Replaces the argument `U≈u(tₙ)` with the next element of the recursive scheme approximating `u(tₙ+dt)` through the formula

```math
u(tₙ+dt)≈e^{-dt D} u(tₙ) + dt \frac{e^{-dt A} - 1}{-dt A} g( u(tₙ) )
```

The matrix `D` should be *diagonal* and the vector of its diagonal values provided together with the nonlinear function `g` by `model`. 
"""
function step!(solver :: EulerExp,
                model :: AbstractModel,
                U  ,
                dt )


    [u1 .= u for (u1,u) in zip(solver.U1,U)]
    model.g!( solver.U1 )
    [u.*=exp.(dt*d) for (d,u) in zip(model.D,U)]
    [u1.*=solver.φ.(dt*d) for (d,u1) in zip(model.D,solver.U1)]
    [u .+= dt * u1 for (u,u1) in zip(U,solver.U1)]

end

"""
    EulerExp_naive()

Runge-Kutta fourth order solver.

A naive version of `Euler`, without argument since no pre-allocation is performed.

"""
struct EulerExp_naive <: TimeSolver
    φ     :: Function
    label :: String

    function EulerExp_naive() 
        φ(z)=(exp(z+eps())-1)/(z+eps())
        new(φ, "exponential Euler (naive)") 
    end
end

"""
    step!(solver :: EulerExp_naive, model :: AbstractModel , U, dt)

A naive version of `step!(solver :: EulerExp, model :: AbstractModel , U, dt)` that does not make use of pre-allocation.
"""
function step!(solver  :: EulerExp_naive,
               model :: AbstractModel ,
               U  ,
               dt )


    U0 = deepcopy(U)
    model.g!( U0 )
    [u.*=exp.(dt*d) for (d,u) in zip(model.D,U)]
    [u0.*=solver.φ.(dt*d) for (d,u0) in zip(model.D,U0)]
    [u .+= dt * u0 for (u,u0) in zip(U,U0)]

end
