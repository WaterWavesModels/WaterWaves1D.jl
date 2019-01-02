export RK4
export step!

"""
    RK4(params)

Runge-Kutta fourth order solver.

"""
mutable struct RK4 <: TimeSolver

    hhat :: Vector{Complex{Float64}}
    uhat :: Vector{Complex{Float64}}
    dh   :: Vector{Complex{Float64}}
    du   :: Vector{Complex{Float64}}

    function RK4( param::NamedTuple )

	n = param.N

        hhat = zeros(Complex{Float64}, n)
        uhat = zeros(Complex{Float64}, n)
        du   = zeros(Complex{Float64}, n)
        dh   = zeros(Complex{Float64}, n)

        new( hhat, uhat, du, dh)

    end

end

function step!(s  :: RK4,
               f  :: AbstractModel,
               h  :: Vector{Complex{Float64}},
               u  :: Vector{Complex{Float64}},
               dt :: Float64)

    s.hhat .= h
    s.uhat .= u
    f( s.hhat, s.uhat)
    s.dh .= s.hhat
    s.du .= s.uhat

    s.hhat .= h .+ dt/2*s.hhat
    s.uhat .= u .+ dt/2*s.uhat
    f( s.hhat, s.uhat)
    s.dh .+= 2 * s.hhat
    s.du .+= 2 * s.uhat

    s.hhat .= h .+ dt/2*s.hhat
    s.uhat .= u .+ dt/2*s.uhat
    f( s.hhat, s.uhat)
    s.dh .+= 2 * s.hhat
    s.du .+= 2 * s.uhat

    s.hhat .= h .+ dt*s.hhat
    s.uhat .= u .+ dt*s.uhat
    f( s.hhat, s.uhat)
    s.dh .+= s.hhat
    s.du .+= s.uhat

    h .+= dt/6 * s.dh
    u .+= dt/6 * s.du
end
