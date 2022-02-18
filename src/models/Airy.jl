export Airy

"""
	Airy(param;label)

Define an object of type `AbstractModel` in view of solving the initial-value problem for
the linear ([Airy](https://en.wikipedia.org/wiki/Airy_wave_theory)) water waves equations.

# Arguments
`param` is of type `NamedTuple` and must contain
- the shallowness parameter `μ` (set the infinite-layer case if `μ=Inf`);
- optionally, `ν` the shallow/deep water scaling factor. By default, `ν=1` if `μ≦1` and `ν=1/√μ` otherwise.
- numerical parameters to construct the mesh of collocation points as `mesh = Mesh(param)`

# Return values
Generate necessary ingredients for solving an initial-value problem via `solve!`:
1. a function `Airy.f!` to be called in explicit time-integration solvers;
2. a function `Airy.mapto` which from `(η,v)` of type `InitialData` provides the raw data matrix on which computations are to be executed;
3. a function `Airy.mapfro` which from such data matrix returns the Tuple of real vectors `(η,v)`, where
    - `η` is the surface deformation;
    - `v` is a velocity variable which is *not* the derivative of the trace of the velocity potential (if not null).

"""
mutable struct Airy <: AbstractModel

	label   :: String
	f!		:: Function
	mapto	:: Function
	mapfro	:: Function

    function Airy(param::NamedTuple; # param is a NamedTuple containing all necessary parameters
		label = "linear (Airy)"  # using a keyword argument allows the user to supersede the default label.
		 )

		# Set up
		μ = param.μ
		if !in(:ν,keys(param)) # set default ν if it is not provided
			ν = min(1,1/√μ)
		else
			ν = param.ν
		end
		# Collocation points and Fourier modes
		m = Mesh(param)
		x, k = m.x, m.k
		# Fourier multipliers
    	∂ₓ	 = 1im * k            # Differentiation
		if μ == Inf || ν == 0
			∂ₓF₁ = 1im * sign.(k)
		else
			∂ₓF₁ = 1im * tanh.(√μ*k)/(√μ*ν)
		end
		# Pre-allocation
		fftη = zeros(Complex{Float64}, m.N)
        fftv = zeros(Complex{Float64}, m.N)

		# Evolution equations are ∂t U = f(U)
		function f!(U)

			fftη .= U[:,1]
			fftv .= U[:,2]

		   	U[:,1] .= -∂ₓF₁.*fftv
		   	U[:,2] .= -∂ₓ.*fftη

		end

		# Build raw data from physical data (discrete Fourier transform)
		function mapto(data::InitialData)
			U = [fft(data.η(x)) fft(data.v(x))]
		end

		# Return physical data `(η,v)` from raw data
		function mapfro(U)
			real(ifft(U[:,1])),real(ifft(U[:,2]))
		end

        new(label, f!, mapto, mapfro )
    end
end
