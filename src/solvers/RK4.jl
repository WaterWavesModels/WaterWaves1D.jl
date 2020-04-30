export RK4
export step!
struct GenericModel <: AbstractModel end
# Useful to have a default keyword argument of type AbstractModel

"""
    RK4(param,model;datasize=2)

Runge-Kutta fourth order solver.

Constructs an object of type `TimeSolver` to be used in `Problem(model, initial, param; solver::TimeSolver)`

`param :: NamedTuple` should contain a value `N` (number of collocation points)

The other arguments are keywords and are optional.

- `datasize :: Int` determines the number of equations solved (default is `datasize=2`);
- `model :: AbstractModel` determines the number of equations solved if it contains a field `:datasize`.

"""
struct RK4 <: TimeSolver

    Uhat :: Array{Complex{Float64},2}
    dU   :: Array{Complex{Float64},2}

    function RK4( param::NamedTuple; model=GenericModel()::AbstractModel, datasize=2::Int )
        if :datasize in fieldnames(typeof(model))
            datasize = model.datasize
        end
        Uhat = zeros(Complex{Float64}, (param.N,datasize))
        dU   = zeros(Complex{Float64}, (param.N,datasize))

        new( Uhat, dU)

    end

end

function step!(s  :: RK4,
               f! :: AbstractModel,
               U  :: Array{Complex{Float64},2},
               dt :: Float64)


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
