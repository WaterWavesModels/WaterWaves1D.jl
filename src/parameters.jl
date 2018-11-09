export Parameters

struct Parameters  

    ϵ  :: Float64
    N  :: Int64
    L  :: Float64
    T  :: Float64
    dt :: Float64

    function Parameters(;kwargs...)

	ϵ  = kwargs[:ϵ]
	N  = kwargs[:N]
	L  = kwargs[:L]
	T  = kwargs[:T]
	dt = kwargs[:dt]
	new( ϵ, N, L, T, dt )

    end

end

