export RK4
export step!

"""
    RK4(arguments;realdata)

Runge-Kutta fourth order solver.

Constructs an object of type `TimeSolver` to be used in `Problem(model, initial, param; solver::TimeSolver)`

Arguments can be either
1. an object of type `AbstractModel`;
2. a `Tuple` `(N,m)` where `N` is the nuber of collocation points and `m` the number of data (equations solved);
3. two integers `N` an `m` as above (the latter is optional, by default `m=2`).
4. a `NamedTuple` containing a key `N` and an integer `m` (the latter is optional, by default `m=2`).

The keyword argument `realdata` is optional, and determines whether pre-allocated vectors are real- or complex-valued.
By default, they are either determined by the model in case `1.`, complex-valued otherwise.

"""
struct RK4 <: TimeSolver

    Uhat :: Array
    dU   :: Array

    function RK4( model :: AbstractModel; realdata=nothing )
        Uhat = mapto(model,Init(x->0*x,x->0*x))
        dU   = copy(Uhat)
        if realdata==true
            Uhat = real.(Uhat);dU = real.(dU)
        end
        if realdata==false
            Uhat = complex.(Uhat);dU = complex.(dU)
        end
        new( Uhat, dU)
    end
    function RK4( datasize; realdata=false )
        if realdata == true
            Uhat = zeros(Float64, datasize)
            dU   = zeros(Float64, datasize)
        else
            Uhat = zeros(Complex{Float64}, datasize)
            dU   = zeros(Complex{Float64}, datasize)
        end
        new( Uhat, dU)
    end
    function RK4( N::Int, datasize=2::Int; realdata=false )
        if realdata == true
            Uhat = zeros(Float64, (N,datasize))
            dU   = zeros(Float64, (N,datasize))
        else
            Uhat = zeros(Complex{Float64}, (N,datasize))
            dU   = zeros(Complex{Float64}, (N,datasize))
        end
        new( Uhat, dU)
    end
    function RK4( param::NamedTuple, datasize=2::Int; realdata=false )
        if realdata == true
            Uhat = zeros(Float64, (param.N,datasize))
            dU   = zeros(Float64, (param.N,datasize))
        else
            Uhat = zeros(Complex{Float64}, (param.N,datasize))
            dU   = zeros(Complex{Float64}, (param.N,datasize))
        end
        new( Uhat, dU)
    end
end

function step!(s  :: RK4,
               f! ,
               U  ,
               dt )


    s.Uhat .= U
    f!( s.Uhat )
    s.dU .= s.Uhat

    s.Uhat .= U .+ dt/2 .* s.Uhat
    f!( s.Uhat )
    s.dU .+= 2 .* s.Uhat

    s.Uhat .= U .+ dt/2 .* s.Uhat
    f!( s.Uhat )
    s.dU .+= 2 .* s.Uhat

    s.Uhat .= U .+ dt .* s.Uhat
    f!( s.Uhat )
    s.dU .+= s.Uhat

    U .+= dt/6 .* s.dU

end
