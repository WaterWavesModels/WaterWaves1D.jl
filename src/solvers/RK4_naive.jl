export RK4_naive
export step!

"""
    RK4_naive(params)

Runge-Kutta fourth order solver.

"""
mutable struct RK4_naive <: TimeSolver

    Uhat :: Array{Complex{Float64},2}
    dU   :: Array{Complex{Float64},2}

    function RK4_naive( param::NamedTuple, model::AbstractModel )

        n = param.N

        Uhat = zeros(Complex{Float64}, (n,model.datasize))
        dU   = zeros(Complex{Float64}, (n,model.datasize))

        new( Uhat, dU)

    end

    function RK4_naive( param::NamedTuple; k=2::Int )

        n = param.N

        Uhat = zeros(Complex{Float64}, (n,k))
        dU   = zeros(Complex{Float64}, (n,k))

        new( Uhat, dU)

    end

end

function step!(s  :: RK4_naive,
               f! :: AbstractModel,
               U  :: Array{Complex{Float64},2},
               dt :: Float64)


    U0 = U
    f!( U0 )
    U1 = U0

    U0 .= U .+ dt/2 .* U1
    f!( U0 )
    U2 = U0

    U0 .= U .+ dt/2 .* U2
    f!( U0 )
    U3 = U0

    U0 .= U .+ dt .* U3
    f!( U0 )
    U4 = U0

    U .+= dt/6 .* (U1 + 2*U2 + 2*U3 + U4 )

end
