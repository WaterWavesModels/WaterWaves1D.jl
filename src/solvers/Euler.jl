export Euler, Euler_naive
export step!

"""
    Euler(arguments;realdata)

Explicit Euler solver.

Construct an object of type `TimeSolver` to be used in `Problem(model, initial, param; solver::TimeSolver)`

Arguments can be either
0. an object of type `AbstractModel`;
1. an `Array` of size `(N,m)` where `N` is the number of collocation points and `m` the number of data (equations solved);
2. a `Tuple` `(N,m)` as above;
3. an integer `N` and an integer `m` as above (the latter is optional, by default `m=2`).
4. a `NamedTuple` containing a key `N` and an integer `m` (the latter is optional, by default `m=2`).

The keyword argument `realdata` is optional, and determines whether pre-allocated vectors are real- or complex-valued.
By default, they are either determined by the model or the type of the array in case `0.` and `1.`, complex-valued otherwise.

"""
struct Euler <: TimeSolver

    U1 :: Array

    function Euler( U :: Array; realdata=nothing )
        U1 = copy(U)
        if realdata==true
            U1 = real.(U1);
        end
        if realdata==false
            U1 = complex.(U1);
        end
        new( U1 )
    end

    function Euler( model :: AbstractModel; realdata=nothing )
        U=model.mapto(Init(x->0*x,x->0*x))
        Euler(U; realdata=realdata)
    end
    function Euler( datasize; realdata=false )
        U = zeros(Float64, datasize)
        Euler(U; realdata=realdata)
    end

    function Euler( N::Int, m=2::Int; realdata=false )
        datasize = (N,m)
        Euler(datasize; realdata=realdata)
    end
    function Euler( param::NamedTuple, datasize=2::Int; realdata=false )
        datasize = (param.N,datasize)
        Euler(datasize; realdata=realdata)
    end
end

function step!(solver :: Euler,
                model :: AbstractModel,
                U  ,
                dt )


    solver.U1 .= U
    model.f!( solver.U1 )
    U .+= dt .* solver.U1

end

"""
    Euler_naive()

Runge-Kutta fourth order solver.

A naive version of `Euler`, without argument since no pre-allocation is performed.

"""
struct Euler_naive <: TimeSolver end

function step!(s  :: Euler_naive,
               model :: AbstractModel ,
               U  ,
               dt )


    U0 = copy(U)
    model.f!( U0 )
    U .+= dt * U0

end
