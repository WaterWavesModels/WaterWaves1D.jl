export NonHydrostatic

"""
    NonHydrostatic(param;kwargs)

Define an object of type `AbstractModel` in view of solving the initial-value problem for
the "Non-hydrostatic" model proposed by [Bristeau, Mangeney, Sainte-Marie and Seguin](https://doi.org/10.3934/dcdsb.2015.20.961)

# Argument
`param` is of type `NamedTuple` and must contain
- dimensionless parameters `ϵ` (nonlinearity) and `μ` (dispersion);
- numerical parameters to construct the mesh of collocation points, if `mesh` is not provided as a keyword argument.

## Optional keyword arguments
- `mesh`: the mesh of collocation points. By default, `mesh = Mesh(param)`;
- `iterative`: solve the elliptic problem through GMRES if `true`, LU decomposition if `false` (default is `true`);
- `precond`: use a (left) preconditioner for GMRES if `true` (default), choose `precond` as the preconditioner if provided;
- `gtol`: relative tolerance of the GMRES algorithm (default is `1e-14`);
- `restart`: the corresponding option of the GMRES algorithm (default is `100`);
- `maxiter`: the corresponding option of GMRES (default is `nothing`);
- `ktol`: tolerance of the Krasny filter (default is `0`, i.e. no filtering);
- `dealias`: dealiasing with Orszag rule `1-dealias/(dealias+2)` (default is `0`, i.e. no dealiasing);
- `label`: a label for future references (default is `"non-hydrostatic"`);

# Return values
Generate necessary ingredients for solving an initial-value problem via `solve!`:
1. a function `NonHydrostatic.f!` to be called in explicit time-integration solvers;
2. a function `NonHydrostatic.mapto` which from `(η,v)` of type `InitialData` provides the raw data matrix on which computations are to be executed;
3. a function `NonHydrostatic.mapfro` which from such data matrix returns the Tuple of real vectors `(η,v,x)`, where
    - `η` is the values of surface deformation at collocation points `x`;
    - `v` is the derivative of the trace of the velocity potential at `x`;
4. additionally, a handy function `NonHydrostatic.mapfrofull` which from data matrix returns the Tuple of real vectors `(η,v,u)`, where
    - `u` corresponds to the layer-averaged velocity.

"""
mutable struct NonHydrostatic <: AbstractModel

	label   :: String
	f!		:: Function
	mapto	:: Function
	mapfro	:: Function
	mapfrofull	:: Function
	info	:: String

    function NonHydrostatic(param::NamedTuple;
							mesh = Mesh(param),
							dealias = 0,
							ktol	= 0,
							iterate = true,
							gtol	= 1e-14,
							precond = true,
							restart = nothing,
							maxiter	= nothing,
							label	= "non-hydrostatic"
							)

		# Set up
		μ 	= param.μ
		ϵ 	= param.ϵ

		if isnothing(maxiter) maxiter = mesh.N end
		if isnothing(restart) restart = min(20,mesh.N) end

		# Print information
		info = "Non-hydrostatic model of Bristeau, Mangeney, Sainte-Marie and Seguin.\n"
		info *= "├─Shallowness parameter μ=$μ, nonlinearity parameter ϵ=$ϵ.\n"
		if dealias == 0
			info *= "├─No dealiasing. "
		else
			info *= "├─Dealiasing with Orszag's rule adapted to power $(dealias + 1) nonlinearity. "
		end
		if ktol == 0
			info *= "No Krasny filter. "
		else
			info *= "Krasny filter with tolerance $ktol."
		end
		if iterate == true
			if precond == false out="out" else out="" end
			info *= "\n└─Elliptic problem solved with GMRES method with$out preconditioning, \
			tolerance $gtol, maximal number of iterations $maxiter, restart after $restart iterations \
			(consider `iterate=false` for non-iterative method). "
		else
			info *= "\n└─Elliptic problem solved with standard LU factorization \
			(consider `iterate=true` for faster results). "
		end
		info *= "\nDiscretized with $(mesh.N) collocation points on [$(mesh.xmin), $(mesh.xmax)]."

		# Pre-allocate useful data
		k = mesh.k
		x 	= mesh.x
		x₀ = mesh.x[1]

		∂ₓ	=  1im * k
		F₁ = 1 ./(1 .+ μ/4*k.^2)
	    if precond == true
			Precond = Diagonal( 1 ./  F₁ )
		elseif precond == false
			Precond = Diagonal( ones(size(k)) )
		else
			Precond = precond
		end
		if dealias == 0
			Π⅔ 	= ones(size(k)) # no dealiasing (Π⅔=Id)
		else
			K = (mesh.kmax-mesh.kmin)/(2+dealias)
			Π⅔ 	= abs.(k) .<= K # Dealiasing low-pass filter
		end
		FFT = exp.(-1im*k*(x.-x₀)');
        IFFT = exp.(1im*k*(x.-x₀)')/length(x);
		M₀ = IFFT * Diagonal( ∂ₓ .* Π⅔)* FFT
		IFFT∂ₓ = IFFT * Diagonal( ∂ₓ .* Π⅔)
        Id = Diagonal(ones(size(x)));
		h = zeros(Complex{Float64}, mesh.N)
		m, u, fftv, fftη, fftu, w = (similar(h),).*ones(6)
		L = similar(FFT)


		# Evolution equations are ∂t U = f(U)
		function f!(U)
			fftη .= U[1]
			h .= 1 .+ ϵ*ifft(fftη)
			fftv .= U[2]
			if iterate == false
				L .= Id - μ/4 * Diagonal(Π⅔) * FFT * Diagonal( 1 ./h ) * M₀ * Diagonal( h.^3 ) * IFFT∂ₓ
				fftu .= L \ fftv
			elseif iterate == true
		        function LL(hatu)
		            hatu- μ/4 *Π⅔.*fft( 1 ./h .* ifft( ∂ₓ .* Π⅔.*fft( h.^3 .* ifft( ∂ₓ .* Π⅔.* hatu ) ) ) )
				end
				fftu .= gmres( LinearMap(LL, length(h); issymmetric=false, ismutating=false) , fftv ;
						restart = restart, maxiter = maxiter, Pl = Precond, reltol = gtol )
			end
			u .= ifft(fftu)
			w .= -h .* ifft(Π⅔.*∂ₓ.*fftu)/2
			m .= ifft(Π⅔.*∂ₓ.*fft(h.*u))
		   	U[1] .= -∂ₓ.*Π⅔.*(fftu .+ ϵ * fft(ifft(fftη) .* u))
			U[2] .= -∂ₓ.*Π⅔.*(fftη .+ ϵ * fft( 1/2*u.^2 ) ) .-
							# The following two lines are equivalent with the last two lines (up to numerical errors)
							#ϵ *μ/2*Π⅔.* fft( 1 ./ h .* ifft(Π⅔.*∂ₓ.*fft(h.*m.*w .+ h.* ifft(Π⅔.*∂ₓ.*fft( h.*u.*w)) ))) .+
							#ϵ *μ/2*Π⅔.* fft( m ./ (h.^2) .* ifft(Π⅔.*∂ₓ.*fft(h.^2 .*w)) )
							#
							ϵ *μ/2*Π⅔.* fft( ϵ *ifft(∂ₓ.*fftη).*u.*ifft(∂ₓ.*fft(w)) .+ w.*ifft(∂ₓ.*∂ₓ.*fft(h.*u)) ) .-
							ϵ *μ/2*Π⅔.* ∂ₓ.*∂ₓ.* fft( h .* u .* w )
			for u in U u[ abs.(u).< ktol ].=0 end
		end

		# Build raw data from physical data.
		# Discrete Fourier transform with, possibly, dealiasing and Krasny filter.
		function mapto(data::InitialData)
			U = [Π⅔ .* fft(data.η(x)), Π⅔ .*fft(data.v(x))]
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
		# Returns `(η,v,u)`, where
		# - `η` is the surface deformation;
		# - `v` is the derivative of the trace of the velocity potential;
		# - `u` corresponds to the layer-averaged velocity.
		# Inverse Fourier transform and take the real part, plus solves the costly elliptic problem for `u`.
		function mapfrofull(U)
				fftη .= U[1]
			   	h .= 1 .+ ϵ*ifft(fftη)
				L .= Id - μ/4 * Diagonal(Π⅔) * FFT * Diagonal( 1 ./h ) * M₀ * Diagonal( h.^3 ) * IFFT∂ₓ

				   real(ifft(U[1])),real(ifft(U[2])),real(ifft(L \ U[2]))
		end

        new(label, f!, mapto, mapfro, mapfrofull, info)
    end
end
