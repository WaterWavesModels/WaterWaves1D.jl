export WhithamGreenNaghdi

"""
    WhithamGreenNaghdi(param;kwargs)

Define an object of type `AbstractModel` in view of solving the initial-value problem for
the fully dispersive Green-Naghdi model proposed by [Duchêne, Israwi and Talhouk](https://doi.org/10.1137/130947064).

# Argument
`param` is of type `NamedTuple` and must contain
- dimensionless parameters `ϵ` (nonlinearity) and `μ` (dispersion);
- numerical parameters to construct the mesh of collocation points as `mesh = Mesh(param)`.

## Optional keyword arguments
- `SGN`: if `true` (default is `false`), compute the Serre-Green-Naghdi (SGN) instead of Whitham-Green-Naghdi (WGN) system (see `SerreGreenNaghdi(param;kwargs)`);
- `iterative`: solve the elliptic problem through GMRES if `true`, LU decomposition if `false` (default is `true`);
- `precond`: use a (left) preconditioner for GMRES if `true` (default), choose `precond` as the preconditioner if provided;
- `gtol`: relative tolerance of the GMRES algorithm (default is `1e-14`);
- `restart`: the corresponding option of the GMRES algorithm (default is `100`);
- `maxiter`: the corresponding option of GMRES (default is `nothing`);
- `ktol`: tolerance of the Krasny filter (default is `0`, i.e. no filtering);
- `dealias`: dealiasing with Orlicz rule `1-dealias/(dealias+2)` (default is `0`, i.e. no dealiasing);
- `label`: a label for future references (default is `"Whitham-Green-Naghdi"`);

# Return values
Generate necessary ingredients for solving an initial-value problem via `solve!`:
1. a function `WhithamGreenNaghdi.f!` to be called in explicit time-integration solvers;
2. a function `WhithamGreenNaghdi.mapto` which from `(η,v)` of type `InitialData` provides the raw data matrix on which computations are to be executed;
3. a function `WhithamGreenNaghdi.mapfro` which from such data matrix returns the Tuple of real vectors `(η,v)`, where
    - `η` is the surface deformation;
    - `v` is the derivative of the trace of the velocity potential;
4. additionally, a handy function `WhithamGreenNaghdi.mapfrofull` which from data matrix returns the Tuple of real vectors `(η,v,u)`, where
    - `u` corresponds to the layer-averaged velocity.

"""
mutable struct WhithamGreenNaghdi <: AbstractModel

	label   :: String
	f!		:: Function
	mapto	:: Function
	mapfro	:: Function
	mapfrofull	:: Function
	info    :: String

    function WhithamGreenNaghdi(param::NamedTuple;
								SGN		= false,
								dealias	= 0,
								ktol	= 0,
								iterate	= true,
								gtol	= 1e-14,
								precond	= true,
								restart	= nothing,
								maxiter	= nothing,
								label	= nothing
								)
		# Set up
		μ 	= param.μ
		ϵ 	= param.ϵ
		mesh = Mesh(param)
		if maxiter == nothing maxiter = mesh.N end
		if restart == nothing restart = min(20,mesh.N) end
		if label == nothing
			if SGN == true
				label = "Serre-Green-Naghdi"
			else
				label = "Whitham-Green-Naghdi"
			end
		end


		# Print information
		info = "$label model.\n"
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
        if SGN == true
			F₁ = 1 ./(1 .+ μ/3*k.^2)
			F₀ = sqrt(μ)*∂ₓ
	    else
			F₁ 	= tanh.(sqrt(μ)*abs.(k))./(sqrt(μ)*abs.(k))
			F₁[1] 	= 1
			F₀ = 1im * sqrt.(3*(1 ./F₁ .- 1)).*sign.(k)
		end
		if precond == true
			#Precond = Diagonal( 1 .+ μ/3*k.^2 )
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
		M₀ = IFFT * Diagonal( F₀ .* Π⅔)* FFT
		IFFTF₀ = IFFT * Diagonal( F₀ .* Π⅔)
        Id = Diagonal(ones(size(x)));
		h = zeros(Complex{Float64}, mesh.N)
		u, fftv, fftη, fftu, hdu = (similar(h),).*ones(5)
		L = similar(FFT)


		# Evolution equations are ∂t U = f(U)
		function f!(U)
			fftη .= U[:,1]
			h .= 1 .+ ϵ*ifft(fftη)
			fftv .= U[:,2]
			if iterate == false
				L .= Id - 1/3 * Diagonal(Π⅔) * FFT * Diagonal( 1 ./h ) * M₀ * Diagonal( h.^3 ) * IFFTF₀
				fftu .= L \ fftv
			elseif iterate == true
		        function LL(hatu)
		            hatu- 1/3 *Π⅔.*fft( 1 ./h .* ifft( F₀ .* Π⅔.*fft( h.^3 .* ifft( F₀ .* Π⅔.* hatu ) ) ) )
				end
				fftu .= gmres( LinearMap(LL, length(h); issymmetric=false, ismutating=false) , fftv ;
						restart = restart, maxiter = maxiter, Pl = Precond, reltol = gtol )
			end
			u .= ifft(fftu)
			hdu .= h .* ifft(Π⅔.*F₀.*fftu)
		   	U[:,1] .= -∂ₓ.*Π⅔.*(fftu .+ ϵ * fft(ifft(fftη) .* u))
			U[:,2] .= -∂ₓ.*Π⅔.*(fftη .+ ϵ * fft( u.*ifft(fftv)
							.- 1/2 * u.^2 .- 1/2 * hdu.^2 ) )
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
		# Return `(η,v,u)`, where
		# - `η` is the surface deformation;
		# - `v` is the derivative of the trace of the velocity potential;
		# - `u` corresponds to the layer-averaged velocity.
		# Inverse Fourier transform and take the real part, plus solves the costly elliptic problem for `u`.
		function mapfrofull(U)
				fftη .= U[:,1]
			   	h .= 1 .+ ϵ*ifft(fftη)
				L .= Id - 1/3 * Diagonal(Π⅔) * FFT * Diagonal( 1 ./h ) * M₀ * Diagonal( h.^3 ) * IFFTF₀

				   real(ifft(U[:,1])),real(ifft(U[:,2])),real(ifft(L \ U[:,2]))
		end

        new(label, f!, mapto, mapfro, mapfrofull, info )
    end
end
