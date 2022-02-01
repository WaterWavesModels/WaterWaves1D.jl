export RK4,RK4_naive
export step!

"""
    RK4(arguments;realdata)

Explicit Runge-Kutta fourth order solver.

Construct an object of type `TimeSolver` to be used in `Problem(model, initial, param; solver::TimeSolver)`

Arguments can be either
0. an object of type `AbstractModel`;
1. an `Array` of size `(N,m)` where `N` is the number of collocation points and `m` the number of data (equations solved);
2. a `Tuple` `(N,m)` as above;
3. an integer `N` and an integer `m` as above (the latter is optional, by default `m=2`).
4. a `NamedTuple` containing a key `N` and an integer `m` (the latter is optional, by default `m=2`).

The keyword argument `realdata` is optional, and determines whether pre-allocated vectors are real- or complex-valued.
By default, they are either determined by the model in case `1.`, complex-valued otherwise.

"""
struct RK4 <: TimeSolver

    U1 :: Array
    dU :: Array
    label :: String

    function RK4( U :: Array; realdata=nothing )
        U1 = copy(U)
        dU = copy(U)
        if realdata==true
            U1 = real.(U1);dU = real.(dU)
        end
        if realdata==false
            U1 = complex.(U1);dU = complex.(dU)
        end
        new( U1, dU, "RK4")
    end

    function RK4( model :: AbstractModel; realdata=nothing )
        U=model.mapto(Init(x->0*x,x->0*x))
        RK4(U; realdata=realdata)
    end
    function RK4( datasize; realdata=false )
        U = zeros(Float64, datasize)
        RK4(U; realdata=realdata)
    end

    function RK4( N::Int, m=2::Int; realdata=false )
        datasize = (N,m)
        RK4(datasize; realdata=realdata)
    end
    function RK4( param::NamedTuple, datasize=2::Int; realdata=false )
        datasize = (param.N,datasize)
        RK4(datasize; realdata=realdata)
    end
end


function step!(s  :: RK4,
                m :: AbstractModel,
                U  ,
                dt )


    s.U1 .= U
    m.f!( s.U1 )
    s.dU .= s.U1

    s.U1 .= U .+ dt/2 .* s.U1
    m.f!( s.U1 )
    s.dU .+= 2 .* s.U1

    s.U1 .= U .+ dt/2 .* s.U1
    m.f!( s.U1 )
    s.dU .+= 2 .* s.U1

    s.U1 .= U .+ dt .* s.U1
    m.f!( s.U1 )
    s.dU .+= s.U1

    U .+= dt/6 .* s.dU

end

"""
    RK4_naive()

Runge-Kutta fourth order solver.

A naive version of `RK4`, without argument since no pre-allocation is performed.

"""
struct RK4_naive <: TimeSolver

    label :: String
    
    function RK4_naive() new("RK4 (naive)") end
end

function step!(s  :: RK4_naive,
               m :: AbstractModel ,
               U  ,
               dt )


    U0 = copy(U)
    m.f!( U0 )
    U1 = copy(U0)

    U0 .= U .+ dt/2 .* U1
    m.f!( U0 )
    U2 = copy(U0)

    U0 .= U .+ dt/2 .* U2
    m.f!( U0 )
    U3 = copy(U0)

    U0 .= U .+ dt .* U3
    m.f!( U0 )
    U4 = copy(U0)

    U .+= dt/6 .* (U1 + 2*U2 + 2*U3 + U4 )

end
