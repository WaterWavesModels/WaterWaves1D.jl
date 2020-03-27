export RK4_naive
export step!

"""
    `RK4_naive()`

Runge-Kutta fourth order solver.
Constructs an object of type `TimeSolver` to be used in `Problem(model, initial, param; solver::TimeSolver)

"""
mutable struct RK4_naive <: TimeSolver
end

function step!(s  :: RK4_naive,
               f! :: AbstractModel,
               U  :: Array{Complex{Float64},2},
               dt :: Float64)


    U0 = U
    f!( U0 )
    U1 = copy(U0)

    U0 .= U .+ dt/2 .* U1
    f!( U0 )
    U2 = copy(U0)

    U0 .= U .+ dt/2 .* U2
    f!( U0 )
    U3 = copy(U0)

    U0 .= U .+ dt .* U3
    f!( U0 )
    U4 = copy(U0)

    U .+= dt/6 .* (U1 + 2*U2 + 2*U3 + U4 )

end
