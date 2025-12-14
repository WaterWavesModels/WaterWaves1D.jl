export Whitham

"""
    Whitham(param; kwargs...)

Define an object of type `AbstractModel` in view of solving the initial-value problem for
two uncoupled Whitham equations, following [Emerald](@cite HoangNguyen2022).

# Argument
`param` is of type `NamedTuple` and must contain
- dimensionless parameters `ϵ` (nonlinearity) and `μ` (dispersion);
- numerical parameters to construct the mesh of collocation points, if `mesh` is not provided as a keyword argument.

## Optional keyword arguments
- `KdV`: if `true` (default is `false`), compute the standard KdV equations instead (see `KdV(param;kwargs)`);
- `BBM`: if `true` (default is `false`), compute the standard BBM equations instead (see `BBM(param;kwargs)`);
- `improved_initial_data`: if `true` (default), improves the naive (first-order) decomposition into right-going and left-going wave;
- `mesh`: the mesh of collocation points. By default, `mesh = Mesh(param)`;
- `ktol`: tolerance of the low-pass Krasny filter (default is `0`, i.e. no filtering);
- `dealias`: dealiasing with Orszag rule `1-dealias/(dealias+2)` (default is `0`, i.e. no dealiasing);
- `label`: a label for future references (default is `"Whitham"`, `"KdV"` or `"BBM"` depending on the equation solved);


# Return values
Generate necessary ingredients for solving an initial-value problem via `solve!`:
1. a function `Whitham.f!` to be called in explicit time-integration solvers;
2. a function `Whitham.mapto` which from `(η,v)` of type `InitialData` provides the raw data matrix on which computations are to be executed.
3. a function `Whitham.mapfro` which from such data matrix returns the Tuple of real vectors `(η,v,x)`, where
  - `η` is the values of surface deformation at collocation points `x`;
  - `v` is the derivative of the trace of the velocity potential at `x`.

"""
mutable struct Whitham <: AbstractModel

	label   :: String
	f!		:: Function
	mapto	:: Function
	mapfro	:: Function
	info 	:: String

    function Whitham(param::NamedTuple;KdV=false,BBM=false,
								improved_initial_data=true,
								mesh = Mesh(param),
								dealias = 0,
								ktol	= 0,
								label 	= nothing
								)

		# Set up
		μ 	= param.μ
		ϵ 	= param.ϵ
		ϵ₀  = ϵ*improved_initial_data

		if KdV == true
			if isnothing(label) label = "KdV" end
		elseif BBM == true
			if isnothing(label) label = "BBM" end
		else
			if isnothing(label) label = "Whitham" end
		end

		# Print information
		info = "$label model.\n"
		info *= "├─Shallowness parameter μ=$μ, nonlinearity parameter ϵ=$ϵ.\n"
		if dealias == 0
			info *= "└─No dealiasing. "
		else
			info *= "└─Dealiasing with Orszag's rule adapted to power $(dealias + 1) nonlinearity. "
		end
		if ktol == 0
			info *= "No Krasny filter. "
		else
			info *= "Krasny filter with tolerance $ktol."
		end
		info *= "\nDiscretized with $(mesh.N) collocation points on [$(mesh.xmin), $(mesh.xmax)]."

		# Pre-allocate data
		x 	= mesh.x
		k = mesh.k
		∂ₓ	=  1im * k
		∂ₓ⁻	=  1 ./(1im * k)
		∂ₓ⁻[1] = 0
		if KdV==true
			F₁ 	= 1 .-μ/6*abs.(k).^2
			F₂  = ones(size(k))
		elseif BBM==true
			F₁ 	= 1 ./(1 .+μ/6*abs.(k).^2)
			F₂  = copy(F₁)
		else
			F₁ 	= tanh.(sqrt(μ)*abs.(k))./(sqrt(μ)*abs.(k))
			F₁[1] 	= 1
			F₁ = sqrt.(F₁)
			F₂ = ones(size(k))

		end
		
		if dealias == 0
			Π⅔ 	= ones(size(k)) # no dealiasing (Π⅔=Id)
		else
			K = (mesh.kmax-mesh.kmin)/(2+dealias)
			Π⅔ 	= abs.(k) .<= K # Dealiasing low-pass filter
		end
		r = zeros(Float64, mesh.N)
        s = zeros(Float64, mesh.N)
		fftr = zeros(Complex{Float64}, mesh.N)
        ffts = zeros(Complex{Float64}, mesh.N)

		# Evolution equations are ∂t U = f(U)
		function f!(U)

		    ffts .= U[2]
			s .= real(ifft(ffts))
			fftr .= U[1]
		   	r .= real(ifft(fftr))

		   	U[1] .= -∂ₓ.*(F₁.*fftr.+3*ϵ/4*Π⅔.*F₂.*fft(r.^2))
		   	U[2] .= ∂ₓ.*(F₁.*ffts.+3*ϵ/4*Π⅔.*F₂.*fft(s.^2))
			for u in U u[ abs.(u).< ktol ].=0 end

		end

		# Build raw data from physical data.
		# Discrete Fourier transform of the suitable variables 
		# with, possibly, dealiasing and Krasny filter.
		function mapto(data::InitialData)
			fftr .= Π⅔ .* ( fft(data.η(x)) + F₁ .* fft(data.v(x)) )/2
			ffts .= Π⅔ .* ( fft(data.η(x)) - F₁ .* fft(data.v(x)) )/2
			r .= real(ifft(fftr))
			s .= real(ifft(ffts))

			U = [fftr - ϵ₀/4*Π⅔ .*fft(ifft(∂ₓ.*fftr).*ifft(∂ₓ⁻.*ffts) + r.*s + 1/2*s.^2) ,
				ffts - ϵ₀/4*Π⅔ .*fft(ifft(∂ₓ.*ffts).*ifft(∂ₓ⁻.*fftr) + r.*s + 1/2*r.^2)  ]
			for u in U u[ abs.(u).< ktol ].=0 end
			return U
		end

		# Reconstruct physical variables from raw data
		# Return `(η,v,x)`, where
		# - `η` is the surface deformation;
		# - `v` is the derivative of the trace of the velocity potential;
		# - `x` is the vector of collocation points
		function mapfro(U)
			ffts .= U[2]
			s .= real(ifft(ffts))
			fftr .= U[1]
		   	r .= real(ifft(fftr))
			newr = fftr + ϵ₀/4*Π⅔ .*fft(ifft(∂ₓ.*fftr).*ifft(∂ₓ⁻.*ffts) + r.*s + 1/2*s.^2)  
			news = ffts + ϵ₀/4*Π⅔ .*fft(ifft(∂ₓ.*ffts).*ifft(∂ₓ⁻.*fftr) + r.*s + 1/2*r.^2)
			real(ifft(newr+news)),real(ifft((newr-news)./F₁)),mesh.x
		end

        new(label, f!, mapto, mapfro, info )
    end
end
