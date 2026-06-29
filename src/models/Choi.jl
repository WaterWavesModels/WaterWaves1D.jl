export Choi

@doc raw"""
    Choi(param; kwargs...)

Define an object of type `AbstractModel` in view of solving the initial-value problem for
(an asymptotically equivalent variant of) the model proposed by [Choi](@cite Choi2022):
```math
  \left\{\begin{array}{l}
∂_tη+∂_x\left(\sum_{m=0}^M  h^{2m+1}\frac{(-μ∂_x^2)^mv }{(2m+1)!}\right) =0,\\[1ex]
\big(	1-\sum_{m=1}^M μ∂_x( h^{2m}∂_x\frac{(-μ∂_x^2)^{m-1}}{(2m)!})\big)∂_tv +∂_xη  \\[1ex] 
		\qquad +∂_x\left(\frac12\sum_{m=0}^M h^{2m}\left(\sum_{j=0}^m \frac{(-μ∂_x^2)^j v}{(2j)!}  \frac{(-μ∂_x^2)^{m-j} v}{(2m-2j)!}-μ\sum_{j=0}^{m-1} \frac{∂_x(-μ∂_x^2)^j v}{(2j+1)!} \frac{∂_x(-μ∂_x^2)^{m-j-1} v}{(2m-2j-1)!} \right)\right),
  \end{array}\right.
```
where ``M`` is the rank of the system, ``h=1 + ϵ η`` , ``η`` the surface deformation and ``v=∂_xϕᵦ`` is the derivative of the trace of the velocity potential at the bottom.

# Argument
`param` is of type `NamedTuple` and must contain
- dimensionless parameters `ϵ` (nonlinearity) and `μ` (dispersion);
- numerical parameters to construct the mesh of collocation points, if `mesh` is not provided as a keyword argument.

## Optional keyword arguments
- `M∈{0,1,2}`: the rank of the system. `M=2` by default;
- `reg`: applies a regularization operator. `reg=false` by default.
- `mesh`: the mesh of collocation points. By default, `mesh = Mesh(param)`;
- `iterative`: solve the elliptic problem through GMRES if `true`, LU decomposition if `false` (default is `true`);
- `precond`: use a (left) preconditioner for GMRES if `true` (default), choose `precond` as the preconditioner if provided;
- `gtol`: relative tolerance of the GMRES algorithm (default is `1e-14`);
- `restart`: the corresponding option of the GMRES algorithm (default is `100`);
- `maxiter`: the corresponding option of GMRES (default is `nothing`);
- `ktol`: tolerance of the Krasny filter (default is `0`, i.e. no filtering);
- `dealias`: dealiasing with Orszag rule `1-dealias/(dealias+2)` (default is `0`, i.e. no dealiasing);
- `label`: a label for future references (default is `"Choi-N"`);

# Return values
Generate necessary ingredients for solving an initial-value problem via `solve!`:
1. a function `Choi.f!` to be called in explicit time-integration solvers;
2. a function `Choi.mapto` which from `(η,v)` of type `InitialData` provides the raw data matrix on which computations are to be executed;
3. a function `Choi.mapfro` which from such data matrix returns the Tuple of real vectors `(η,v,x)`, where
  - `η` is the values of surface deformation at collocation points `x`;
  - `v` is the derivative of the trace of the velocity potential;
4. additionally, a handy function `Choi.mapfrofull` which from data matrix returns the Tuple of real vectors `(η,v,u,vᵦ)`, where
  - `v` is the derivative of the trace of the velocity potential;
  - `u` corresponds to the layer-averaged horizontal velocity;
  - `vᵦ` corresponds to the horizontal velocity at the bottom.

"""
mutable struct Choi <: AbstractModel

	label   :: String
	f!		:: Function
	mapto	:: Function
	mapfro	:: Function
	mapfrofull	:: Function
	info    :: String

    function Choi(param::NamedTuple; M=2, reg = false,
								mesh = Mesh(param),
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

		if isnothing(maxiter) maxiter = mesh.N end
		if isnothing(restart) restart = min(20,mesh.N) end
		if isnothing(label)
			label = "Choi-$M"
		end


		# Print information
		info = "$label model.\n"
		info *= "├─Shallowness parameter μ=$μ, nonlinearity parameter ϵ=$ϵ.\n"
		info *= "├─Order M=$M."
		if reg == true  
			info *= " Regularization added in the model.\n"
		else
			info *= " No regularization added in the model.\n"
		end
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
		Dμ  = -μ*(∂ₓ.^2)
		dμ	= sqrt(μ)*∂ₓ

		if precond == true
			p = ones(size(k))
			for m in 0:M
				p.+= (Dμ.^m)/factorial(m)
			end
			Precond = Diagonal( p )
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
	    Id = Diagonal(ones(size(x)));
		h = zeros(Complex{Float64}, mesh.N)
		q, fftv, fftη, fftu = (similar(h),).*ones(4)
		L = similar(FFT)
		Z = zero(FFT)

		


		function f1(h,fftv)
			y=h.*ifft(fftv)
			for m in 1:M
				y+=h.^(2*m+1).*ifft((Dμ.^m).*fftv)/factorial(2*m+1)
			end
			return fft(y)
		end

		function f2(h,fftv)
			y=zero(fftv)
			for m in 1:M
				y+=h.^(2*m).*ifft(dμ.*(Dμ.^(m-1)).*fftv)/factorial(2*m)
			end
			return fftv-dμ.*fft(y)
		end

		function f3(h,fftv)
			y=ifft(fftv).^2
			z=similar(y)
			for m in 1:M
				z.=ifft((Dμ.^m).*fftv).*ifft(fftv)/factorial(2*m)
				for j in 0:m-1
				z+=ifft((Dμ.^j).*fftv).*ifft((Dμ.^(m-j)).*fftv)/factorial(2*j)/factorial(2*m-2*j)
				z-=ifft(dμ.*(Dμ.^j).*fftv).*ifft(dμ.*(Dμ.^(m-j-1)).*fftv)/factorial(2*j+1)/factorial(2*m-2*j-1)
				end
				y+=h.^(2*m).*z
			end
			return 1/2*fft(y)
		end


		function regul(f) 
			ifft(exp.(-μ*k.^2).*fft(f))
			ifft(min.(1,1 ./(μ*k.^2)).*fft(f))
		end




		# Evolution equations are ∂t U = f(U)
		function f!(U)
			fftη .= U[1]
			h .= 1 .+ ϵ*ifft(fftη)
			fftv .= U[2]

			U[1] .= -∂ₓ.*Π⅔.*(f1(h,fftv))
			q      .= -∂ₓ.*Π⅔.*(fftη .+ ϵ * f3(h,fftv) )

			if reg == true
				h = regul(h)
			end

			if iterate == false
				Z .= zero(FFT)
				for m in 1:M
					Z+=Diagonal(h.^(2*m))*IFFT*Diagonal(dμ.*(Dμ.^(m-1)))/factorial(2*m)
				end
				L = Id-Diagonal(dμ)*FFT*Z
				fftu .= L \ q
			elseif iterate == true
				function F2(fftv) f2(h,fftv) end
				fftu .= gmres( LinearMap(F2, length(h); issymmetric=false, ismutating=false) , q ;
						restart = restart, maxiter = maxiter, Pl = Precond, reltol = gtol )
			end

			U[2] .= fftu
			for u in U u[ abs.(u).< ktol ].=0 end
		end

		# Build raw data from physical data.
		# Discrete Fourier transform with, possibly, dealiasing and Krasny filter.
		function mapto(data::InitialData)

			fftη .= fft(data.η(x))
			h .= 1 .+ ϵ*ifft(fftη)

			if reg == true
				h = regul(h)
			end

			 
			if iterate == false
				Z .= zero(FFT)
				for m in 1:M
					Z+=Diagonal(h.^(2*m))*IFFT*Diagonal(dμ.*(Dμ.^(m-1)))/factorial(2*m)
				end
				L = Id-Diagonal(dμ)*FFT*Z
				fftu .= L \ fft(data.v(x))
			elseif iterate == true
				function F2(fftv) f2(h,fftv) end
				fftu .= gmres( LinearMap(F2, length(h); issymmetric=false, ismutating=false) , fft(data.v(x)) ;
						restart = restart, maxiter = maxiter, Pl = Precond, reltol = gtol )
			end


			U = [Π⅔.*fftη, Π⅔.*fftu]
			for u in U u[ abs.(u).< ktol ].=0 end
			return U
		end

		# Reconstruct physical variables from raw data
		# Return `(η,v,x)`, where
		# - `η` is the surface deformation;
		# - `v` is the derivative of the trace of the velocity potential;
		# - `x` is the vector of collocation points
		function mapfro(U)
			fftη .= U[1]
			h .= 1 .+ ϵ*ifft(fftη)
			fftv .= U[2]

			real(ifft(U[1])),real(ifft(f2(h,fftv))),mesh.x
		end
		# Return `(η,v,u,vb)`, where
		# - `η` is the surface deformation;
		# - `v` is the derivative of the trace of the velocity potential;
		# - `u` corresponds to the layer-averaged velocity.
		# - `vᵦ` corresponds to the horizontal velocity at the bottom.

		# Inverse Fourier transform and take the real part, plus solves the costly elliptic problem for `u`.
		function mapfrofull(U)
				fftη .= U[1]
			   	h .= 1 .+ ϵ*ifft(fftη)
				fftv .= U[2]

				real(ifft(U[1])),real(ifft(f2(h,fftv))),real(ifft(f1(h,fftv))./h),real(ifft(fftv))
		end

        new(label, f!, mapto, mapfro, mapfrofull, info )
    end
end


