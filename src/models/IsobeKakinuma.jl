export IsobeKakinuma

"""
    IsobeKakinuma(param;kwargs)

Define an object of type `AbstractModel` in view of solving the initial-value problem for
the Isobe-Kakinuma model proposed by [Isobe](https://doi.org/10.1061/9780784400890.023).

# Argument
`param` is of type `NamedTuple` (or a collection `NamedTuple`s) of and must contain
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
- `label`: a label for future references (default is `"Isobe-Kakinuma"`);

# Return values
Generate necessary ingredients for solving an initial-value problem via `solve!`:
1. a function `IsobeKakinuma.f!` to be called in explicit time-integration solvers;
2. a function `IsobeKakinuma.mapto` which from `(η,v)` of type `InitialData` provides the raw data matrix on which computations are to be executed;
3. a function `IsobeKakinuma.mapfro` which from such data matrix returns the Tuple of real vectors `(η,v)`, where
    - `η` is the surface deformation;
    - `v` is the derivative of the trace of the velocity potential;
4. additionally, a handy function `IsobeKakinuma.mapfrofull` which from data matrix returns the Tuple of real vectors `(η,v,Φ)`, where
    - `Φ` is the Vector of the basis functions `ϕi` (`i∈{0,...,N}`).

"""
mutable struct IsobeKakinuma <: AbstractModel

	label   :: String
	f!		:: Function
	mapto	:: Function
	mapfro	:: Function
	mapfrofull	:: Function
	info	:: String

    function IsobeKakinuma(param::NamedTuple;
				dealias = 0,
				ktol	= 0,
				iterate = true,
				gtol	= 1e-14,
				precond = true,
				restart	= nothing,
				maxiter	= nothing,
				label	= "Isobe-Kakinuma"
				)

		# Set up
		μ 	= param.μ
		ϵ 	= param.ϵ
		mesh = Mesh(param)

		if isnothing(maxiter) maxiter = mesh.N end
		if isnothing(restart) restart = min(20,mesh.N) end

		# Print information
		info = "Isobe-Kakinuma model of order 2.\n"
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
		if dealias == 0
			Π⅔ 	= ones(size(k)) # No dealiasing (Π⅔=Id)
		else
			K = (mesh.kmax-mesh.kmin)/(2+dealias)
			Π⅔ 	= abs.(k) .<= K # Dealiasing low-pass filter
		end
		FFT = exp.(-1im*k*(x.-x₀)');
        IFFT = exp.(1im*k*(x.-x₀)')/length(x);
        Id = Diagonal(ones(size(x)));

		if precond == true
			Precond = lu([  Id   μ*Id  ;
					1/2*Diagonal( ∂ₓ.^2) (Id + μ/10 * Diagonal(∂ₓ.^2 .* Π⅔) ) ])
		elseif precond == false
			Precond = Diagonal( ones( size( [k;k])) )
		else
			Precond = precond
		end

		h = zeros(Complex{Float64}, mesh.N)
		u, fftv, fftη, u, w, G = (similar(h),).*ones(6)
		C = similar([h ; h])
		guess = 0*C
		fftϕ = reshape(C,:,2)
		L = similar([FFT FFT ; FFT FFT])


		# Evolution equations are ∂t U = f(U)
		function f!(U)
			fftη .= U[:,1]
			h .= 1 .+ ϵ*ifft(fftη)
			fftv .= U[:,2]
			C .= [ zero(fftv) ; -1/2 * ∂ₓ.* fftv]
			if iterate == false
				L .= [Id    μ*(FFT * Diagonal( h.^2 ) * IFFT .* Π⅔)  ;
						1/2*Diagonal( ∂ₓ.^2 .* Π⅔) (Id + μ/10 * Diagonal(Π⅔) * FFT * Diagonal( h.^2 ) * IFFT * Diagonal( ∂ₓ.^2 .* Π⅔)) ]

				fftϕ .= reshape(L \ C,:,2)
			elseif iterate == true # does not work yet
				guess .= fftϕ[:]  # bof, le guess n'ameliore pas tellement les performances
		        function LL(fftϕ)
		            [fftϕ[1:end÷2] + μ*fft( (h.^2) .* ifft(  Π⅔ .* fftϕ[end÷2+1:end]))   ;
					1/2*∂ₓ.^2 .* Π⅔.* fftϕ[1:end÷2] .+ fftϕ[end÷2+1:end] .+ μ/10*Π⅔.*fft( (h.^2) .* ifft( ∂ₓ.^2 .* Π⅔ .* fftϕ[end÷2+1:end] ) ) ]
				end
				fftϕ .= reshape( gmres!( guess, LinearMap(LL, 2*length(x); issymmetric=false, ismutating=false) , C ;
						restart = restart, maxiter = maxiter, Pl = Precond, reltol = gtol ) , : ,2)
			end
			u .= ifft( ∂ₓ.* fftϕ[:,1] .+ fftv) .+ μ * (h.^2)  .* ifft(∂ₓ.* fftϕ[:,2])
			w .= 2*h.* ifft( fftϕ[:,2])
			G .= -∂ₓ.*Π⅔.* fft( h.* ifft( ∂ₓ.* fftϕ[:,1] .+ fftv) .+ μ*(h.^3)/3  .* ifft(∂ₓ.* fftϕ[:,2]) )

		   	U[:,1] .= G
			U[:,2] .= -∂ₓ .* (fftη .+ ϵ * Π⅔ .* fft( μ * w .* ifft(-G)
								.+ 1/2 * (u.^2 .+ μ * w.^2 ) ) )
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
		# Returns `(η,v,u)`, where
		# - `η` is the surface deformation;
		# - `v` is the derivative of the trace of the velocity potential;
		# - `u` corresponds to the layer-averaged velocity.
		# Inverse Fourier transform and take the real part, plus solves the costly elliptic problem for `u`.
		function mapfrofull(U)   # remains to be done
				fftη .= U[:,1]
			   	h .= 1 .+ ϵ*ifft(fftη)
				L .= Id - 1/3 * Diagonal(Π⅔) * FFT * Diagonal( 1 ./h ) * M₀ * Diagonal( h.^3 ) * IFFTF₀

				   real(ifft(U[:,1])),real(ifft(U[:,2])),real(ifft(L \ U[:,2]))
		end

        new(label, f!, mapto, mapfro, mapfrofull, info)
    end
end
