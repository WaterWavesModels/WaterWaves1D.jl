export Parameters

"""
   Parameters( ϵ  = 1/2,
               N  = 2^12,
               L  = 10,
               T  = 5,
               dt = 0.001)

"""
struct Parameters

    N   :: Int64
    L   :: Float64
    T   :: Float64
    dt  :: Float64
	dtr :: Float64
	arg :: NamedTuple

    function Parameters(;kwargs...)

	N   = kwargs[:N]
	L   = kwargs[:L]
	T   = kwargs[:T]
	dt  = kwargs[:dt]
	if :dtr in keys(kwargs)
		dtr = kwargs[:dtr]
	else
		dtr = kwargs[:dt]
	end
	arg = values(kwargs)
	new( N, L, T, dt, dtr, arg )

    end

end
