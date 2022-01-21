export SquareRootDepth

"""
    SquareRootDepth(param;kwargs)

Define an object of type `AbstractModel` in view of solving the initial-value problem for
the "√D" model proposed by [Cotter, Holm and Percival](https://doi.org/10.1098/rspa.2010.0124)

# Argument
`param` is of type `NamedTuple` and must contain
- dimensionless parameters `ϵ` (nonlinearity) and `μ` (dispersion);
- numerical parameters to construct the mesh of collocation points as `mesh = Mesh(param)`.

## Optional keyword arguments
- `iterative`: solve the elliptic problem through GMRES if `true`, LU decomposition if `false` (default is `true`);
- `precond`: use a (left) preconditioner for GMRES if `true` (default), choose `precond` as the preconditioner if provided;
- `gtol`: relative tolerance of the GMRES algorithm (default is `1e-14`);
- `restart`: the corresponding option of the GMRES algorithm (default is `100`);
- `maxiter`: the corresponding option of GMRES (default is `nothing`);
- `ktol`: tolerance of the Krasny filter (default is `0`, i.e. no filtering);
- `dealias`: dealiasing with Orlicz rule `1-dealias/(dealias+2)` (default is `0`, i.e. no dealiasing);
- `label`: a label for future references (default is `"square-root depth"`);
- `verbose`: prints information if `true` (default is `true`).

# Return values
Generate necessary ingredients for solving an initial-value problem via `solve!`:
1. a function `SquareRootDepth.f!` to be called in explicit time-integration solvers;
2. a function `SquareRootDepth.mapto` which from `(η,v)` of type `InitialData` provides the raw data matrix on which computations are to be executed;
3. a function `SquareRootDepth.mapfro` which from such data matrix returns the Tuple of real vectors `(η,v)`, where
    - `η` is the surface deformation;
    - `v` is the derivative of the trace of the velocity potential;
4. additionally, a handy function `SquareRootDepth.mapfrofull` which from data matrix returns the Tuple of real vectors `(η,v,u)`, where
    - `u` corresponds to the layer-averaged velocity.

"""
mutable struct SquareRootDepth <: AbstractModel

	label   :: String
	f!		:: Function
	mapto	:: Function
	mapfro	:: Function
	mapfrofull	:: Function

    function SquareRootDepth(param::NamedTuple;
							dealias = 0,
							ktol	= 0,
							iterate	= true,
							gtol	= 1e-14,
							precond	= true,
							restart	= nothing,
							maxiter	= nothing,
							label	= "square-root depth",
							verbose	= true)

		if verbose @info "Build the √D model of Cotter, Holm and Percival." end

		μ 	= param.μ
		ϵ 	= param.ϵ
		mesh = Mesh(param)

		k = mesh.k
		x 	= mesh.x
		x₀ = mesh.x[1]

		∂ₓ	=  1im * mesh.k
		F₁ = 1 ./(1 .+ μ/3*k.^2)
	    if precond == true
			Precond = Diagonal( 1 ./  F₁ )
		elseif precond == false
			Precond = Diagonal( ones(size(k)) )
		else
			Precond = precond
		end
		K = mesh.kmax * (1-dealias/(2+dealias))
		Π⅔ 	= abs.(mesh.k) .<= K # Dealiasing low-pass filter
		if dealias == 0
			if verbose @info "no dealiasing" end
			Π⅔ 	= ones(size(mesh.k))
		elseif verbose
			@info "dealiasing : spectral scheme for power  $(dealias + 1) nonlinearity"
		end
		if iterate == true && verbose
			@info "elliptic problem solved with GMRES method"
		elseif verbose
			@info "elliptic problem solved with LU decomposition"
		end
		FFT = exp.(-1im*k*(x.-x₀)');
        IFFT = exp.(1im*k*(x.-x₀)')/length(x);
		M₀ = IFFT * Diagonal( ∂ₓ .* Π⅔) * FFT
		IFFT∂ₓFFT = IFFT * Diagonal( ∂ₓ .* Π⅔) * FFT
        Id = Diagonal(ones(size(x)));
		h = zeros(Complex{Float64}, mesh.N)
		m, u, fftv, fftη, fftu, w = (similar(h),).*ones(6)
		L = similar(FFT)
		if maxiter == nothing maxiter = mesh.N end
		if restart == nothing restart = min(20,mesh.N) end


		# Evolution equations are ∂t U = f!(U)
		function f!(U)
			fftη .= U[:,1]
			h .= 1 .+ ϵ*ifft(fftη)
			fftv .= U[:,2]
			if iterate == false
				L .= Id - μ/3 * Diagonal( ∂ₓ .* Π⅔) * FFT * Diagonal( 1 ./h ) * IFFT∂ₓFFT * Diagonal( h ) * IFFT
				fftu .= L \ fftv
			elseif iterate == true
		        function LL(hatu)
		            hatu - μ/3 * ∂ₓ .* Π⅔.*fft( 1 ./h .* ifft( ∂ₓ .* Π⅔ .* fft( h .* ifft(hatu ) ) ) )
				end
				fftu .= gmres( LinearMap(LL, length(h); issymmetric=false, ismutating=false) , fftv ;
						restart = restart, maxiter = maxiter, Pl = Precond, reltol = gtol )
			end
			u .= ifft(fftu)
			w .= ifft(Π⅔.*∂ₓ.*fft(h.*u))
		   	U[:,1] .= -∂ₓ.*Π⅔.*(fftu .+ ϵ * fft(ifft(fftη) .* u))
			U[:,2] .= -∂ₓ.*Π⅔.*(fftη .+ ϵ/2 * fft( u.^2 )
							.+ ϵ*μ/6 * fft( (w./h).^2 ) )
			U[abs.(U).< ktol ].=0
		end

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
		# Returns `(η,v,u)`, where
		# - `η` is the surface deformation;
		# - `v` is the derivative of the trace of the velocity potential;
		# - `u` corresponds to the layer-averaged velocity.
		# Inverse Fourier transform and take the real part, plus solves the costly elliptic problem for `u`.
		function mapfrofull(U)
				fftη .= U[:,1]
			   	h .= 1 .+ ϵ*ifft(fftη)
				L .=  Id - μ/3 * Diagonal( ∂ₓ .* Π⅔) * FFT * Diagonal( 1 ./h ) * IFFT∂ₓFFT * Diagonal( h ) * IFFT

				   real(ifft(U[:,1])),real(ifft(U[:,2])),real(ifft(L \ U[:,2]))
		end

        new(label, f!, mapto, mapfro, mapfrofull)
    end
end
