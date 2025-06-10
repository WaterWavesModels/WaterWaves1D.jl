export RK4,RK4_naive
export step!

"""
    RK4(arguments;realdata)

Explicit Runge-Kutta fourth order solver.

Construct an object of type `TimeSolver` to be used in `Problem(model, initial, param; solver::TimeSolver)`

Arguments can be either
0. an object of type `AbstractModel`, typically the model you will solve with the solver;
1. an `Array` which has the size of the objects that the solver will manipulate (typically a vector of `systemsize` elements of size `N` where `N` is the number of collocation points and `systemsize` the number of solved equations);
2. a `(datasize,systemsize)` where `datasize` is the size of scalar variables (typically `N` the number of collocation points) and `datasize` (optional, by default `systemsize=2`) the number of solved equations);
3. `(param,systemsize)` where `param` is a `NamedTuple` containing a key `N` describing the number of collocation points, and `systemsize` the number of solved equations (optional, by default `systemsize=2`).

The keyword argument `realdata` is optional, and determines whether pre-allocated vectors are real- or complex-valued.
By default, they are either determined by the model or the type of the array in case `0.` and `1.`, complex-valued in case `2.`.

"""
struct RK4 <: TimeSolver

    U1 :: Array
    dU :: Array
    label :: String

    function RK4( U :: Array; realdata=nothing )
        U1 = deepcopy(U)
        dU = deepcopy(U)
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
        RK4( U; realdata=realdata)
    end
    function RK4( param::NamedTuple, systemsize=2::Int; realdata=nothing )
        RK4( [Array{Complex{Float64}}(undef,param.N) for _ in 1:systemsize]  ; realdata=realdata)
    end
    function RK4( datasize, systemsize=2::Int; realdata=nothing )
        RK4( [Array{Complex{Float64}}(undef,datasize) for _ in 1:systemsize] ; realdata=realdata)
    end
end

using RecursiveArrayTools

function step!(s  :: RK4,
    m :: AbstractModel,
    U  ,
    dt )

    u = VectorOfArray(U)
    u1 = VectorOfArray(s.U1)
    du = VectorOfArray(s.dU)

    u1 .= u 

    m.f!( s.U1 )

    du .= u1 

    u1 .= u .+ dt/2 .* u1 

    m.f!( s.U1 )

    du .+= 2 .* u1 

    u1 .= u .+ dt/2 .* u1 

    m.f!( s.U1 )

    du .+= 2 .* u1 

    u1 .= u .+ dt .* u1 

    m.f!( s.U1 )

    du .+= u1 

    u .+= dt/6 .* du 

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


    U0 = deepcopy(U)
    m.f!( U0 )
    U1 = deepcopy(U0)

    [u0 .= u .+ dt/2 .* u1 for (u0,u,u1) in zip(U0,U,U1)]
    m.f!( U0 )
    U2 = deepcopy(U0)

    [u0 .= u .+ dt/2 .* u2 for (u0,u,u2) in zip(U0,U,U2)]
    m.f!( U0 )
    U3 = deepcopy(U0)

    [u0 .= u .+ dt .* u3 for (u0,u,u3) in zip(U0,U,U3)]
    m.f!( U0 )
    U4 = deepcopy(U0)

    [u .+= dt/6 .* (u1 + 2*u2 + 2*u3 + u4 ) for (u,u1,u2,u3,u4) in zip(U,U1,U2,U3,U4)]

end