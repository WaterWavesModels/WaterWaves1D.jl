export WhithamBoussinesq

"""
    WhithamBoussinesq(param;kwargs)

Define an object of type `AbstractModel` in view of solving the initial-value problem for
a Boussinesq-type model with full-dispersion property.

# Argument
`param` is of type `NamedTuple` and must contain
- dimensionless parameters `ϵ` (nonlinearity) and `μ` (dispersion);
- numerical parameters to construct the mesh of collocation points as `mesh = Mesh(param)`


## Optional keyword arguments
- `Boussinesq`: if `true` (default is `false`), compute the standard Boussinesq system instead (see `Boussinesq(param;kwargs)`);
- a parameter `α` which determines the model solved:
    - If `α = 1` (default), then the model has been introduced in [Dinvay, Dutykh and Kalisch](https://doi.org/10.1016/j.apnum.2018.09.016);
    - If `α = 1/2`, then the model is a quasilinear version;
    - If `α < 1/2`, then expect instabilities stemming from ill-posedness of the model.
- `ktol`: tolerance of the low-pass Krasny filter (default is `0`, i.e. no filtering);
- `dealias`: dealiasing with Orlicz rule `1-dealias/(dealias+2)` (default is `0`, i.e. no dealiasing);
- `label`: a label for future references (default is `"Whitham-Boussinesq"`);
- `verbose`: prints information if `true` (default is `true`).


# Return values
Generate necessary ingredients for solving an initial-value problem via `solve!`:
1. a function `WhithamBoussinesq.f!` to be called in explicit time-integration solvers;
2. a function `WhithamBoussinesq.mapto` which from `(η,v)` of type `InitialData` provides the raw data matrix on which computations are to be executed.
3. a function `WhithamBoussinesq.mapfro` which from such data matrix returns the Tuple of real vectors `(η,v)`, where
    - `η` is the surface deformation;
    - `v` is the derivative of the trace of the velocity potential.

"""
mutable struct WhithamBoussinesq <: AbstractModel

	label   :: String
	f!		:: Function
	mapto	:: Function
	mapfro	:: Function

    function WhithamBoussinesq(param::NamedTuple;Boussinesq=false,
								α = 1, a = -1/3, b = 1/3,
								dealias = 0,
								ktol	= 0,
								label 	= nothing,
								verbose	=true)

		# Set up
		μ 	= param.μ
		ϵ 	= param.ϵ

		if Boussinesq == true
			if label == nothing label = "Boussinesq" end
			info_param = "a=$a, b=$b, c=0 and d=$b"
		else
			if label == nothing label = "Whitham-Boussinesq" end
			info_param = "α=$α"
		end

		if verbose # Print information
			info = "Build the $label model with $info_param.\n"
			info *= "Shallowness parameter μ=$μ, nonlinearity parameter ϵ=$ϵ.\n"
			if dealias == 0
				info *= "No dealiasing. "
			else
				info *= "Dealiasing with Orszag's rule adapted to power $(dealias + 1) nonlinearity. "
			end
			if ktol == 0
				info *= "No Krasny filter. "
			else
				info *= "Krasny filter with tolerance $ktol."
			end
			info *= "\nShut me up with keyword argument `verbose = false`."
			@info info
		end

		# Pre-allocate data
		mesh = Mesh(param)
		x 	= mesh.x
		k = mesh.k
		∂ₓ	=  1im * k
		if Boussinesq==false
			F₁ 	= tanh.(sqrt(μ)*abs.(k))./(sqrt(μ)*abs.(k))
			F₁[1] 	= 1
			F₂ = F₁.^α
		else
			F₂ = 1 ./(1 .+μ*b*abs.(k).^2)
			F₁ 	= (1 .-μ*a*abs.(k).^2).*(F₂.^2)
		end

		if dealias == 0
			Π⅔ 	= ones(size(k)) # no dealiasing (Π⅔=Id)
		else
			K = (mesh.kmax-mesh.kmin)/(2+dealias)
			Π⅔ 	= abs.(k) .<= K # Dealiasing low-pass filter
		end
		η = zeros(Float64, mesh.N)
        v = zeros(Float64, mesh.N)
		fftη = zeros(Complex{Float64}, mesh.N)
        fftv = zeros(Complex{Float64}, mesh.N)

		# Evolution equations are ∂t U = f(U)
		function f!(U)

		    fftv .= U[:,2]
			fftη .= F₂.*U[:,2]
			v .= real(ifft(fftη))
			fftη .= U[:,1]
		   	η .= real(ifft(U[:,1]))

		   	U[:,1] .= -∂ₓ.*(F₁.*fftv.+ϵ*Π⅔.*F₂.*fft(η.*v))
		   	U[:,2] .= -∂ₓ.*(fftη.+ϵ/2*Π⅔.*fft(v.^2))
			U[abs.(U).< ktol ].=0

		end

		# Build raw data from physical data.
		# Discrete Fourier transform with, possibly, dealiasing and Krasny filter.
		function mapto(data::InitialData)
			U = [Π⅔ .* fft(data.η(x)) Π⅔ .*fft(data.v(x))]
			U[abs.(U).< ktol ].=0
			return U
		end

		# Return `(η,v)`, where
		# - `η` is the surface deformation;
		# - `v` is the derivative of the trace of the velocity potential.
		# Inverse Fourier transform and takes the real part.
		function mapfro(U)
			real(ifft(U[:,1])),real(ifft(U[:,2]))
		end

        new(label, f!, mapto, mapfro)
    end
end
