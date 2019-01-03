export RK4
export step!

"""
    RK4(params)

Runge-Kutta fourth order solver.

"""
mutable struct RK4 <: TimeSolver

    Uhat :: Array{Complex{Float64},2}
    dU 	 :: Array{Complex{Float64},2}

    function RK4( param::NamedTuple, k::Int )

		n = param.N

        Uhat = zeros(Complex{Float64}, (n,k))
		dU 	 = zeros(Complex{Float64}, (n,k))

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

    s.Uhat .= U .+ dt/2*s.Uhat
    f!( s.Uhat )
    s.dU .+= 2 * s.Uhat

    s.Uhat .= U .+ dt/2*s.Uhat
	f!( s.Uhat )
    s.dU .+= 2 * s.Uhat

    s.Uhat .= U .+ dt*s.Uhat
	f!( s.Uhat )
    s.dU .+= s.Uhat

    U .+= dt/6 * s.dU
end
