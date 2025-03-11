export threescale

"""
    threescale(param;kwargs)

Define an object of type `AbstractModel` in view of solving the initial-value problem for
a three-scale problem 

# Argument
`param` is of type `NamedTuple` and must contain
- dimensionless parameters `ϵ` (nonlinearity) and `δ` (dispersion);
- the relaxation parameter `a`;
- numerical parameters to construct the mesh of collocation points, if `mesh` is not provided as a keyword argument.

## Optional keyword arguments
- `simple`: if `true` (default is `true`), the model is simple, otherwise it is closer to LCT;
- `str`: if `true` (default is `true`), the model has an extra structure, otherwise it does not;
- `prep`: `∈{0,1,2}` and represent the level of preparation of the initial data (default is `1`);
- `mesh`: the mesh of collocation points. By default, `mesh = Mesh(param)`;
- `iterative`: solve the elliptic problem (to construct initial data) through GMRES if `true`, LU decomposition if `false` (default is `true`);
- `precond`: use a (left) preconditioner for GMRES if `true` (default), choose `precond` as the preconditioner if provided;
- `gtol`: relative tolerance of the GMRES algorithm (default is `1e-14`);
- `restart`: the corresponding option of the GMRES algorithm (default is `100`);
- `maxiter`: the corresponding option of GMRES (default is `nothing`);
- `ktol`: tolerance of the Krasny filter (default is `0`, i.e. no filtering);
- `dealias`: dealiasing with Orszag rule `1-dealias/(dealias+2)` (default is `0`, i.e. no dealiasing);
- `label`: a label for future references;

# Return values
Generate necessary ingredients for solving an initial-value problem via `solve!`:
1. a function `threescale.f!` to be called in explicit time-integration solvers;
2. a function `threescale.mapto` which from `(η,v)` of type `InitialData` provides the raw data matrix on which computations are to be executed;
3. a function `threescale.mapfro` which from such data matrix returns the Tuple of real vectors `(η,v,x)`, where
	- `η` is the values of surface deformation at collocation points `x`;
	- `v` is the derivative of the trace of the velocity potential at `x`;
4. additionally, a handy function `threescale.mapfrofull` which from data matrix returns the Tuple of real vectors `(η,v,u,p,w)`, where
	- `u` corresponds to the layer-averaged horizontal velocity.
	- `p` corresponds to the relaxed (artificial) layer-averaged non-hydrostatic pressure;
	- `w` corresponds to the relaxed (artificial) layer-averaged vertical velocity.

"""
mutable struct threescale <: AbstractModel

	label   :: String
	f!		:: Function
	mapto	:: Function
	mapfro	:: Function
	mapfrofull	:: Function
	info    :: String

    function threescale(param::NamedTuple; 
								simple = [1,1,1,1],
								str = true,
								prep = 1,
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
		a	= param.a
		δ 	= param.δ
		μ   = δ^2
		ϵ 	= param.ϵ

		if isnothing(maxiter) maxiter = mesh.N end
		if isnothing(restart) restart = min(20,mesh.N) end
		if isnothing(label)
			if simple == true name = "simple" else name = "" end
			if str == true
				label = "$name 3-scale model with structure"
			else
				label = "$name 3-scale without structure"
			end
		end


		# Print information
		info = "$label model.\n"
		info *= "├─Relaxation parameter a=$a.\n"
		info *= "├─Shallowness parameter δ=$δ, nonlinearity parameter ϵ=$ϵ.\n"
		info *= "├─Initial data prepared of order $prep.\n"
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

		if precond == true
			Precond = Diagonal( 1 .+ μ/3*k.^2 )
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
		Dx = IFFT * Diagonal( Π⅔ .* ∂ₓ)* FFT
        Id = Diagonal(ones(size(x)));
		h = zeros(Complex{Float64}, mesh.N)
		η, u, p, w, fftη, fftu, fftv, fftp, fftw = (similar(h),).*ones(9)
		L = similar(FFT)


		# Evolution equations are ∂t U = f(U)
		function f!(U;a=a)
			if simple != false
				fftη .= U[:,1]; fftu .= U[:,2];
				fftp .= U[:,3]; fftw .= U[:,4];
				η .= ifft(fftη)
				if str == false
					η .+= ifft(fftu)/2
				end
	
				U[:,1] .= -simple[1]*∂ₓ.* fftu
							#-∂ₓ.*Π⅔.*fft(h .* u)
				U[:,2] .=  - Π⅔.*fft(1 ./(1 .+ ϵ*simple[2]*η) .* ifft(∂ₓ.* (simple[1]*fftη .+ a*δ *fftp) ) )
							#-∂ₓ.*Π⅔.*(fftη .+ ϵ/2 * fft( u.^2) ) - a*μ * Π⅔.*fft(1 ./h .* ifft(∂ₓ.* fft( hFG.*p ) ) )
				U[:,3] .= -a*Π⅔.*fft(1 ./(1 .+ ϵ*simple[3]*η) .*ifft(2*fftw.+ δ * ∂ₓ.*fftu) )
							#  -a*Π⅔.*fft((2*w.+hFG.*ifft(∂ₓ.*fftu))./h) - ϵ*Π⅔.*fft( u .* ifft(∂ₓ.* fftp)  )
				U[:,4] .= 3/2*a*Π⅔.*fft(1 ./(1 .+ ϵ*simple[4]*η) .*ifft(fftp)) #- ϵ*Π⅔.*fft( u .* ifft(∂ₓ.* fftw ) )
				U[abs.(U).< ktol ].=0
			else
				fftη .= U[:,1]
				h .= 1 .+ ϵ*ifft(fftη)
				fftu .= U[:,2]; u.= ifft(fftu);
				fftp .= U[:,3]; p.= ifft(fftp);
				fftw .= U[:,4]; w.= ifft(fftw);
				if str == false
					h .+= ϵ*ifft(fftu)/2
				end

				U[:,1] .= -∂ₓ.*fftu
							#-∂ₓ.*Π⅔.*fft(h .* u)
				U[:,2] .= -∂ₓ.*fftη - a*δ * Π⅔.*fft(1 ./h .* ifft(∂ₓ.* fft( h.*p ) ) )
							#-∂ₓ.*Π⅔.*(fftη .+ ϵ/2 * fft( u.^2) ) - a*μ * Π⅔.*fft(1 ./h .* ifft(∂ₓ.* fft( hFG.*p ) ) )
				U[:,3] .= -a*Π⅔.*fft((2*w.+ δ * h.*ifft(∂ₓ.*fftu))./h) 
							#  -a*Π⅔.*fft((2*w.+hFG.*ifft(∂ₓ.*fftu))./h) - ϵ*Π⅔.*fft( u .* ifft(∂ₓ.* fftp)  )
				U[:,4] .= 3/2*a*Π⅔.*fft(p./h) #- ϵ*Π⅔.*fft( u .* ifft(∂ₓ.* fftw ) )
				U[abs.(U).< ktol ].=0
			end
		end

		# Build raw data from physical data.
		# Discrete Fourier transform with, possibly, dealiasing and Krasny filter.
		function mapto(data::InitialData)
			fftη .= fft(data.η(x)) 
			fftv .= fft(data.v(x)) 
			fftu .= copy(fftv)

			h .= 1 .+ ϵ*ifft(fftη)
			if str == false
				η .+= ifft(fftu)/2
			end

			# if iterate == false
			# 	L .= Id - μ/3 * Diagonal(Π⅔) * FFT * Diagonal( 1 ./h ) * Dx * Diagonal( h.^3 ) * IFFT * Diagonal( Π⅔ .* ∂ₓ )
			# 	fftu .= L \ fftv
			# else#if iterate == true
		    #     function LL(hatu)
		    #         hatu- μ/3 *Π⅔.*fft( 1 ./h .* ifft( Π⅔ .* ∂ₓ .*fft( h.^3 .* ifft( Π⅔ .* ∂ₓ .* hatu ) ) ) )
			# 	end
			# 	fftu .= gmres( LinearMap(LL, length(h); issymmetric=false, ismutating=false) , fftv ;
			# 			restart = restart, maxiter = maxiter, Pl = Precond, reltol = gtol )
			# end
			if prep >= 2
				# if iterate == false
				# 	L .= FFT* Diagonal( 3 ./h.^3 ) * IFFT - μ * Diagonal(  Π⅔ .* ∂ₓ ) * FFT *  Diagonal( 1 ./h ) * IFFT * Diagonal( Π⅔ .* ∂ₓ )
				# 	fftp .= fft(ifft(Π⅔ .*(L \ (Π⅔ .* ( ∂ₓ.^2 .*fftη + 2 * fft(ifft(Π⅔ .* ∂ₓ.* fftu).^2)))) / a)./h)
				# elseif iterate == true
				# 	function ll(hatu)
				# 		Π⅔.*fft( 3 ./h.^3 .* ifft(hatu))- μ * Π⅔ .* ∂ₓ .*fft( 1 ./h .* ifft( Π⅔ .*∂ₓ .* hatu ) ) 
				# 	end
				# 	fftp .= fft(ifft(Π⅔ .* gmres( LinearMap(ll, length(h); issymmetric=false, ismutating=false) , Π⅔ .* (  ∂ₓ.^2 .*fftη + 2 * fft(ifft(Π⅔ .* ∂ₓ.* fftu).^2)) ;
				# 			restart = restart, maxiter = maxiter, Pl = Precond, reltol = gtol ) / a) ./h)
				# elseif iterate == 1/2
					U₀ = [Π⅔ .* fftη Π⅔ .* fftu fft(zero(x)) Π⅔ .* (-1/2*fft(h.*ifft(∂ₓ.*fftu)))]
					dt = param.dt.^2
					# Calcul U(dt)
					U0 = copy(U₀)
    				f!( U0 ; a = a^2)
    				U1 = copy(U0)

    				U0 .= U₀ .+ dt/2 .* U1
    				f!( U0 ; a = a^2)
    				U2 = copy(U0)

    				U0 .= U₀ .+ dt/2 .* U2
    				f!( U0 ; a = a^2)
    				U3 = copy(U0)

    				U0 .= U₀ .+ dt .* U3
    				f!( U0 ; a = a^2)
    				U4 = copy(U0)

    				U₊ = U₀ + dt/6 .* (U1 + 2*U2 + 2*U3 + U4 )

					# Calcul U(-dt)
					U0 = copy(U₀)
    				f!( U0 ; a = a^2)
    				U1 = -copy(U0)

    				U0 .= U₀ .+ dt/2 .* U1
    				f!( U0 ; a = a^2)
    				U2 = -copy(U0)

    				U0 .= U₀ .+ dt/2 .* U2
    				f!( U0 ; a = a^2)
    				U3 = -copy(U0)

    				U0 .= U₀ .+ dt .* U3
    				f!( U0 ; a = a^2)
    				U4 = -copy(U0)

    				U₋ = U₀ + dt/6 .* (U1 + 2*U2 + 2*U3 + U4 )

					d2th = ifft(U₊[:,1]+U₋[:,1]-2*U₀[:,1])/dt^2
					dth = ifft(U₊[:,1]-U₋[:,1])/(2*dt)
					dtdxh = ifft(∂ₓ.*U₊[:,1]-∂ₓ.*U₋[:,1])/(2*dt)
					d2xh = ifft(∂ₓ.*∂ₓ.*U₀[:,1])
					dtu = ifft(U₊[:,2]-U₋[:,2])/(2*dt)
					dxu = ifft(∂ₓ.*U₀[:,2])
					dxh = ifft(∂ₓ.*U₀[:,1])
					u = ifft(U₀[:,2])

					# Need to relaod everything since this has been modified by f!
					fftη .= fft(data.η(x)) 
					fftv .= fft(data.v(x)) 
					h .= 1 .+ ϵ*ifft(fftη)
					fftu .= fft(u)

					fftp = δ/3/a*Π⅔ .*fft(  h .* (d2th+dtu.*dxh+2*u.*dtdxh+u.*dxu.*dxh+u.*u.*d2xh) )


				#end
				
			else
				fftp.=fft(zero(x))
			end
			if prep == -1 
				fftw .= +1/2*∂ₓ.*fftu # fft(zero(x))
			elseif prep == 0 
					fftw .=  fft(zero(x))
			else 
				if simple != false
					fftw .=-1/2*∂ₓ.*fftu
				else
					fftw .=-1/2*fft(h.*ifft(∂ₓ.*fftu))
				end
			end
			# U = [Π⅔ .* fftη Π⅔ .* fftu Π⅔ .* fftp Π⅔ .* fftw]
			U = [Π⅔ .* fftη Π⅔ .* fftu Π⅔ .* (δ*fftp) Π⅔ .* (δ*fftw)]

			U[abs.(U).< ktol ].=0
			return U
		end

		# Reconstruct physical variables from raw data
		# Return `(η,v,x)`, where
		# - `η` is the surface deformation;
		# - `v` is the derivative of the trace of the velocity potential;
		# - `x` is the vector of collocation points
		function mapfro(U)
			fftη .= U[:,1]
			#h .= 1 .+ ϵ*ifft(fftη)
			fftu .= U[:,2]
			fftv .= fftu #- μ/3 *Π⅔.*fft( 1 ./h .* ifft(  Π⅔.* ∂ₓ .* fft( h.^3 .* ifft(  Π⅔.* ∂ₓ .* fftu ) ) ) )
			fftp .= U[:,3]; 
			fftw .= U[:,4];
			real(ifft(fftp)),real(ifft(fftv)),mesh.x
		end
		# Return `(η,v,u,p,w,x)`, where
		# - `η` is the surface deformation;
		# - `v` is the derivative of the trace of the velocity potential;
		# - `u` corresponds to the layer-averaged horizontal velocity.
		# - `p` corresponds to the relaxed (artificial) layer-averaged non-hydrostatic pressure;
		# - `w` corresponds to the relaxed (artificial) layer-averaged vertical velocity.
		# - `x` is the vector of collocation points
		function mapfrofull(U)
			fftη .= U[:,1]
			#h .= 1 .+ ϵ*ifft(fftη)
			fftu .= U[:,2]
			fftv .= fftu #- μ/3 *Π⅔.*fft( 1 ./h .* ifft(  Π⅔.* ∂ₓ .* fft( h.^3 .* ifft(  Π⅔.* ∂ₓ .* fftu ) ) ) )
			real(ifft(fftη)),real(ifft(fftv)),real(ifft(fftu)),real(ifft(U[:,3])),real(ifft(U[:,4])),mesh.x
		end

        new(label, f!, mapto, mapfro, mapfrofull, info )
    end


end
