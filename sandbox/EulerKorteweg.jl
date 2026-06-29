export EulerKorteweg,EulerKorteweg_GP

@doc raw"""
    EulerKorteweg(param; kwargs...)

Define an object of type `AbstractModel` in view of solving the initial-value problem for the
Euler-Korteweg model
```math
  \left\{\begin{array}{l}
  ∂_tη+∂_x(hv)=0,\\[1ex]
  ∂_tv+∂_xη+ϵv∂_xv=\frac{κ^2}{2ϵ}∂_x(\frac{∂_x^2(√h)}{√h}),
  \end{array}\right.
```
where ``h = 1+ϵη``.

# Argument
`param` is of type `NamedTuple` and must contain
- the dimensionless parameter `ϵ` (nonlinearity);
- the dimensionless parameter `κ` (capillarity);
- numerical parameters to construct the mesh of collocation points, if `mesh` is not provided as a keyword argument.

## Optional keyword arguments
- `mesh`: the mesh of collocation points. By default, `mesh = Mesh(param)`;
- `ktol`: tolerance of the low-pass Krasny filter (default is `0`, i.e. no filtering);
- `dealias`: no dealisasing if set to `0` or `false` (default), otherwise `1/(3*dealias)` modes are set to `0` (corresponding to standard 2/3 Orszag rule if `dealias` is set to `1` or `true`);
- `smooth`: A smooth low-pass filter (whose scaling is defined by ) if set to `0` or `false` (default), otherwise only `2/(3*dealias)*(1-smooth/2)` modes are kept untouched;
- `label`: a label for future references (default is `"Saint-Venant"`);

# Return values
Generate necessary ingredients for solving an initial-value problem via `solve!`:
1. a function `EulerKorteweg.f!` to be called in explicit time-integration solvers;
2. a function `EulerKorteweg.mapto` which from `(η,v)` of type `InitialData` provides the raw data matrix on which computations are to be executed;
3. a function `EulerKorteweg.mapfro` which from such data matrix returns the Tuple of real vectors `(η,v,x)`, where
    - `η` is the values of surface deformation at collocation points `x`;
    - `v` is the layer-averaged velocity (or the derivative of the trace of the velocity potential) at `x`.

See also [`EulerKorteweg_GP`](@ref) for the method consisting in solving the Gross-Pitaevskii equation, 
and [`EulerKorteweg_Grenier`](@ref) for the formulation introduced in Grenier.
"""
mutable struct EulerKorteweg <: AbstractModel

	label 	:: String
	f!		:: Function
	mapto	:: Function
	mapfro	:: Function
	info	:: String

    function EulerKorteweg(param::NamedTuple;
						mesh = Mesh(param),
						dealias=0,smooth=false,
						ktol=0,
						label="Euler-Korteweg"
						)

		# Set up
		ϵ 	= param.ϵ
		κ 	= param.κ

		info = "Euler-Korteweg model.\n"
		info *= "├─Nonlinearity parameter ϵ=$(param.ϵ).\n"
		if dealias == 0
			info *= "└─No dealiasing. "
		elseif dealias == 1
			info *= "├─Dealiasing with Orszag's rule adapted to power 2 nonlinearity. \n"
		else
			info *= "├─Dealiasing at user-defined mode. \n"
		end
		if smooth != 0
			info *= "└─Lipschitz low-pass filter is applied. "
		end
		if ktol == 0
			info *= "No Krasny filter. "
		else
			info *= "Krasny filter with tolerance $ktol."
		end
		info *= "\nDiscretized with $(mesh.N) collocation points on [$(mesh.xmin), $(mesh.xmax)]."

		# Pre-allocate useful data
		k  = mesh.k
		x  = mesh.x
		∂ₓ	=  1im * k

		if dealias == 0
			Π⅔ 	= ones(size(k)) # no dealiasing (Π⅔=Id)
			Π 	= ones(size(k)) # no dealiasing (Π⅔=Id)
		else
			K = (mesh.kmax-mesh.kmin)/(3*dealias)
			Π⅔ 	= abs.(k) .<= K # Dealiasing low-pass filter
			if smooth ==0
				Π 	= abs.(k) .<= K # Dealiasing low-pass filter
			else
				Π  = max.( 0 , min.( 1 , 2/smooth*( 1 .-abs.(k)/K) )).^2	
			end
		end
		η = zeros(Float64, mesh.N)
        v = zeros(Float64, mesh.N)
		fftη = zeros(Complex{Float64}, mesh.N)
        fftv = zeros(Complex{Float64}, mesh.N)
		sqrth = zeros(Float64, mesh.N)
		

		# Evolution equations are ∂t U = f(U)
		function f!(U)
			fftη .= U[1]
			fftv .= U[2]
			η .= real(ifft(fftη))
			v .= real(ifft(fftv))
			sqrth .= sqrt.(abs.(1 .+ ϵ*η))

			U[1] .= -∂ₓ.*( fftv .+ ϵ * Π.* fft( η .* v) )
			U[2] .= -∂ₓ.*( fftη .+ ϵ * Π.* fft( v.^2 )/2 .- 
						κ^2/ϵ*fft(ifft((∂ₓ.^2).*fft(sqrth))./sqrth)/2 )
			for u in U u[ abs.(u).< ktol ].=0 end
		end

		# Build raw data from physical data.
		# Discrete Fourier transform with, possibly, dealiasing and Krasny filter.
		function mapto(data::InitialData)
			U = [Π⅔.* fft(data.η(x)), Π⅔.*fft(data.v(x))]
			for u in U u[ abs.(u).< ktol ].=0 end
			return U
		end

		# Reconstruct physical variables from raw data
		# Return `(η,v,x)`, where
		# - `η` is the surface deformation;
		# - `v` is the derivative of the trace of the velocity potential;
		# - `x` is the vector of collocation points
		function mapfro(U)
			real(ifft(U[1])),real(ifft(U[2])),mesh.x
		end

		new(label, f!, mapto, mapfro, info)
    end
end

"""
	EulerKorteweg_GP(param;kwargs)

Same as [`EulerKorteweg`](@ref), but using the Madelung transform to formulate the Euler-Korteweg equation as the Gross-Pitaevskii equation:
``ψ = √h e^{i ϵϕ/κ}`` where ``∂ₓϕ = v`` satisfies ``iκ∂ₜψ + κ^2 ∂_x^2 ψ = (|ψ|^2-1)ψ =0``.
"""
mutable struct EulerKorteweg_GP <: AbstractModel

	label 	:: String
	f!		:: Function
	D		:: AbstractArray
	g!		:: Function
	mapto	:: Function
	mapfro	:: Function
	info	:: String

    function EulerKorteweg_GP(param::NamedTuple;
						mesh = Mesh(param),
						dealias=0,smooth=false,
						ktol=0,
						label="Euler-Korteweg"
						)

		# Set up
		ϵ 	= param.ϵ
		κ 	= param.κ

		info = "Saint-Venant model.\n"
		info *= "├─Nonlinearity parameter ϵ=$(param.ϵ).\n"
		if dealias == 0
			info *= "└─No dealiasing. "
		elseif dealias == 1
			info *= "├─Dealiasing with Orszag's rule adapted to power 2 nonlinearity. \n"
		else
			info *= "├─Dealiasing at user-defined mode. \n"
		end
		if smooth != 0
			info *= "└─Lipschitz low-pass filter is applied. "
		end
		if ktol == 0
			info *= "No Krasny filter. "
		else
			@error("The Krasny filter will *not* be applied. Use 'EulerKorteweg' instead")
		end
		info *= "\nDiscretized with $(mesh.N) collocation points on [$(mesh.xmin), $(mesh.xmax)]."

		# Pre-allocate useful data
		k  = mesh.k
		x  = mesh.x
		∂ₓ	=  1im * k

		if dealias == 0
			Π⅔ 	= ones(size(k)) # no dealiasing (Π⅔=Id)
			Π 	= ones(size(k)) # no dealiasing (Π⅔=Id)
		else
			K = (mesh.kmax-mesh.kmin)/(3*dealias)
			Π⅔ 	= abs.(k) .<= K # Dealiasing low-pass filter
			if smooth ==0
				Π 	= abs.(k) .<= K # Dealiasing low-pass filter
			else
				Π  = max.( 0 , min.( 1 , 2/smooth*( 1 .-abs.(k)/K) )).^2	
			end
		end
		ψ = zeros(Complex{Float64}, mesh.N)
		fftψ = zeros(Complex{Float64}, mesh.N)

		# Evolution equations are ∂t U = f(U)
		function f!(U)
			fftψ .= U[1]
			ψ .= ifft(fftψ)

			U[1] .= -1im/κ*Π.* fft(ψ.*(abs.(ψ).^2 .-1) ) + 1im*κ/2 *(∂ₓ.^2).*fftψ
			for u in U u[ abs.(u).< ktol ].=0 end
		end

		# Evolution equations rewritten as ∂t U = Diag(D) U + g(U)
		D = [1im*κ/2 *(∂ₓ.^2),]
		
		function g!(U)
			fftψ .= U[1]
			ψ .= ifft(fftψ)

			U[1] .= -1im/κ*Π.* fft(ψ.*(abs.(ψ).^2 .-1) ) 
			for u in U u[ abs.(u).< ktol ].=0 end
		end

		# Build raw data from physical data.
		# Apply Madelung transform, then
		# discrete Fourier transform with, possibly, dealiasing and Krasny filter.
		function mapto(data::InitialData)
			h = 1 .+ ϵ*data.η(x)
			fftϕ = fft(data.v(x))./∂ₓ
			fftϕ[1] = complex(0)
			ψ=sqrt.(abs.(h)).*exp.(1im*ϵ/κ*real.(ifft(fftϕ)))
			U = Π⅔ .* fft(ψ)
			U[ abs.(U).< ktol ].=0 
			return [U,]
		end

		# Reconstruct physical variables from raw data
		# Return `(η,v,x)`, where
		# - `η` is the surface deformation;
		# - `v` is the derivative of the trace of the velocity potential;
		# - `x` is the vector of collocation points
		function mapfro(U)
			fftψ=U[1]
			(abs.(ifft(fftψ)).^2 .-1)/ϵ,κ/ϵ*	imag((ifft(∂ₓ.*fftψ)).*conj(ifft(fftψ)))./(abs.(ifft(fftψ)).^2),
			#κ*real.(ifft(∂ₓ.*fft(angle.(ifft(U))))),
			mesh.x
		end

		new(label, f!, D, g!, mapto, mapfro, info)
    end
end

"""
	EulerKorteweg_Grenier(param;kwargs)

Same as [`EulerKorteweg`](@ref), but using the system used by Grenier in the context of the semiclassical limit of the Gross-Pitaevskii equation.

"""
mutable struct EulerKorteweg_Grenier <: AbstractModel

	label 	:: String
	f!		:: Function
	mapto	:: Function
	mapfro	:: Function
	info	:: String

    function EulerKorteweg_Grenier(param::NamedTuple;
						mesh = Mesh(param),
						dealias=0,smooth=false,
						ktol=0,
						label="Euler-Korteweg"
						)

		# Set up
		ϵ 	= param.ϵ
		κ 	= param.κ

		info = "Euler-Korteweg model.\n"
		info *= "├─Nonlinearity parameter ϵ=$(param.ϵ).\n"
		if dealias == 0
			info *= "└─No dealiasing. "
		elseif dealias == 1
			info *= "├─Dealiasing with Orszag's rule adapted to power 2 nonlinearity. \n"
		else
			info *= "├─Dealiasing at user-defined mode. \n"
		end
		if smooth != 0
			info *= "└─Lipschitz low-pass filter is applied. "
		end
		if ktol == 0
			info *= "No Krasny filter. "
		else
			info *= "Krasny filter with tolerance $ktol."
		end
		info *= "\nDiscretized with $(mesh.N) collocation points on [$(mesh.xmin), $(mesh.xmax)]."

		# Pre-allocate useful data
		k  = mesh.k
		x  = mesh.x
		∂ₓ	=  1im * k

		if dealias == 0
			Π⅔ 	= ones(size(k)) # no dealiasing (Π⅔=Id)
			Π 	= ones(size(k)) # no dealiasing (Π⅔=Id)
		else
			K = (mesh.kmax-mesh.kmin)/(3*dealias)
			Π⅔ 	= abs.(k) .<= K # Dealiasing low-pass filter
			if smooth ==0
				Π 	= abs.(k) .<= K # Dealiasing low-pass filter
			else
				Π  = max.( 0 , min.( 1 , 2/smooth*( 1 .-abs.(k)/K) )).^2	
			end
		end
		v = zeros(Complex{Float64}, mesh.N)
		a = zeros(Complex{Float64}, mesh.N)
		ffta = zeros(Complex{Float64}, mesh.N)
		fftv = zeros(Complex{Float64}, mesh.N)
		

		# Evolution equations are ∂t U = f(U)
		function f!(U)
			ffta .= U[1]
			fftv .= U[2]
			a .= ifft(ffta)
			v .= real(ifft(fftv))

			U[1] .= 1im*κ/2* (∂ₓ.^2).* ffta - ϵ * Π.* fft( a.*ifft(∂ₓ.*fftv)/2 .+ v.*ifft(∂ₓ.*ffta) )
			U[2] .= -∂ₓ.*(  ϵ * Π.* fft( v.^2 )/2 .+ Π.* fft( abs.(a).^2 .-1 )/ϵ )
			for u in U u[ abs.(u).< ktol ].=0 end
		end

		# Build raw data from physical data.
		# discrete Fourier transform with, possibly, dealiasing and Krasny filter.
		function mapto(data::InitialData)
			h = 1 .+ ϵ*data.η(x)
			U = [Π⅔.* fft(sqrt.(abs.(h))), Π⅔.*fft(data.v(x))]
			for u in U u[ abs.(u).< ktol ].=0 end
			return U
		end
		# Reconstruct physical variables from raw data
		# Return `(η,v,x)`, where
		# - `η` is the surface deformation;
		# - `v` is the derivative of the trace of the velocity potential;
		# - `x` is the vector of collocation points
		function mapfro(U)
			fftϕ = U[2]./∂ₓ
			fftϕ[1] = complex(0)
			ϕ=real.(ifft(fftϕ))
			fftψ = fft(ifft(U[1]).*exp.(1im/κ*ϵ*ϕ))
			(abs.(ifft(fftψ)).^2 .-1)/ϵ,κ/ϵ*	imag((ifft(∂ₓ.*fftψ)).*conj(ifft(fftψ)))./(abs.(ifft(fftψ)).^2),mesh.x
		end

		new(label, f!, mapto, mapfro, info)
    end
end
