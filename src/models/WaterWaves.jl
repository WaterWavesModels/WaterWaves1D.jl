export WaterWaves

"""
    WaterWaves(param; kwargs)

Define an object of type `AbstractModel` in view of solving the initial-value problem for
the water waves system (via conformal mapping, see [Zakharov, Dyachenko and Vasilyev](https://doi.org/10.1016/S0997-7546(02)01189-5)).

# Argument
`param` is of type `NamedTuple` and must contain
- dimensionless parameters `ϵ` (nonlinearity) and `μ` (dispersion);
- optionally, `ν` the shallow/deep water scaling factor. By default, `ν=1` if `μ≦1` and `ν=1/√μ` otherwise. Set the infinite-layer case if `ν=0`, or `μ=Inf`.
- numerical parameters to construct the mesh of collocation points as `mesh = Mesh(param)`.

## Optional keyword arguments
- `IL`: Set the infinite-layer case if `IL=true`, in which case `ϵ` is the steepness parameter. Default is `false`.
- `method ∈ {1,2,3}`: method used to initialize the conformal mapping, as a fix-point problem `F(u)=u`
    - if `method == 1`, use standard contraction fix-point iteration;
    - if `method == 2`, use Newton algorithm with GMRES iterative solver to invert the Jacobian;
    - if `method == 3`, use Newton algorithm with direct solver to invert the Jacobian;
- `tol`: (relative) tolerance of the fix-point algorithm (default is `1e-16`);
- `maxiter`: the maximal number of iteration in the fix-point algorithm (default is `100`);
- `ktol`: tolerance of the low-pass Krasny filter (default is `0`, i.e. no filtering);
- `dealias`: dealiasing with Orlicz rule `1-dealias/(dealias+2)` (default is `0`, i.e. no dealiasing);
- `label`: a label for future references (default is `"water waves"`);
- `verbose`: prints information if `true` (default is `true`).

# Return values
Generate necessary ingredients for solving an initial-value problem via `solve!`:
1. a function `WaterWaves.f!` to be called in the explicit time-integration solver (also `WaterWaves.f1!` and `WaterWaves.f2!` for the symplectic Euler solver);
2. a function `WaterWaves.mapto` which from `(η,v)` of type `InitialData` provides the raw data matrix on which computations are to be executed;
3. a function `WaterWaves.mapfro` which from such data matrix returns the Tuple of real vectors `(x,η,v)`, where
    - `x` is a vector of collocation points (non-regularly spaced);
    - `η` is the surface deformation at points `x`;
    - `v` is the derivative of the trace of the velocity potential at points `x`.

"""
mutable struct WaterWaves <: AbstractModel

	label   :: String
	f!		:: Function
	f1!		:: Function
	f2!		:: Function
	mapto	:: Function
	mapfro	:: Function
	info	:: String

    function WaterWaves(param::NamedTuple;
						IL	    = false,
						method  = 1,
						tol	    = 1e-16,
						maxiter = 100,
						dealias	= 0,
						ktol	= 0,
						label	= "water waves",
						verbose	= true)

		# Set up
		μ 	= param.μ
		ϵ 	= param.ϵ
		if !in(:ν,keys(param))
			if μ > 1
				ν = 1/sqrt(μ)
				nu = "1/√μ (deep water case)"
			else
				ν = 1
				nu = "1 (shallow water case)"
			end
		else
			ν = param.ν
			nu = "$ν"
		end
		if μ == Inf || ν==0 # infinite layer case
			IL = true;  # IL (=Infinite layer) is a flag to be used in functions tanh,cotanh,xdcotanh
			μ = 1; ν = 1; # Then we should set μ=ν=1 in subsequent formula.
		end
		mesh = Mesh(param)


		# Print information
		info = "Water waves system.\n"
		if IL == true
			info *= "├─Steepness parameter ϵ=$ϵ (infinite depth case).\n"
		else
			info *= "├─Shallowness parameter μ=$μ, nonlinearity parameter ϵ=$ϵ, \
					scaling parameter ν=$nu.\n"
		end
		info *= "├─Initial data built with method $method: "
		if method == 1 info *= "standard contraction fix-point iteration." end
		if method == 2 info *= "Newton algorithm with GMRES iterative solver to invert the Jacobian." end
		if method == 3 info *= "Newton algorithm with direct solver to invert the Jacobian." end
		info *= " Relative tolerance $tol and maximum $maxiter iterations.\n"
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

		# Pre-allocate useful data

		x 	= mesh.x
		x₀ = x[1]
		k 	= mesh.k

    	∂ₓ	=  1im * k             # Differentiation Fourier multiplier
		if dealias == 0
			Π⅔ 	= ones(size(k)) # no dealiasing (Π⅔=Id)
		else
			K = (mesh.kmax-mesh.kmin)/(2+dealias)
			Π⅔ 	= abs.(k) .<= K # Dealiasing low-pass filter
		end

		fz = zeros(Complex{Float64}, mesh.N)
        z = zeros(Float64, mesh.N)
        phi = zeros(Float64, mesh.N)
		ξ = zeros(Float64, mesh.N)
		xv = zeros(Float64, mesh.N)
		Dxv = zeros(Float64, mesh.N)
		Dz = zeros(Float64, mesh.N)
		Dphi = zeros(Float64, mesh.N)
		J = zeros(Float64, mesh.N)
		M1 = zeros(Float64, mesh.N)
		M2 = zeros(Float64, mesh.N)
		q0 = 0.

		# some useful functions
		if IL == true
			cotanh   = x -> sign.(x)
			xdcotanh = x -> zero(x)
			mytanh   = x -> sign.(x)
		else
			cotanh   = x -> (x.!=0) ./ tanh.(x+(x.==0))
			xdcotanh = x -> -x ./ (sinh.(x+(x.==0)).^2)
			mytanh 	 = x -> tanh.(x)
		end

		function meanf(f) #computes the mean of a function, being given its FFT
			real(f[1])/length(f)
		end
		function mean(f) #computes the mean of a function
			sum(f)/length(f)
		end

		# Build raw data from physical data.
		# Use a conformal change of coordinate
		function mapto(data::InitialData)
		  if data.η(x)==zero(x) && data.v(x)==zero(x)  #useful when defining the type of initial data in some solvers
			  return zeros(Complex{Float64}, (mesh.N,2))
		  else
			# preallocate some variables to save memory use
			fη = zeros(Complex{Float64}, mesh.N)
			u0 = zeros(Complex{Float64}, mesh.N)
			nu0 = zeros(Complex{Float64}, mesh.N)
			z0 = zeros(Float64, mesh.N)
			v0 = zeros(Float64, mesh.N)
			dx = zeros(Float64, mesh.N)

			# We will solve F(u)=0
			function F( u )
				dx.=real.(ifft(u))
				fη .= fft(data.η(x+ϵ/ν*dx))
			    return u+ν*sqrt(μ)*1im*Π⅔ .*cotanh(sqrt(μ)*(1+ϵ*meanf(fη))*k).*fη
			end

			# The Newton alogorithm demands the Jacobian of F,
			# which involves the derivative of η
			hatdη=∂ₓ.*fft(data.η(x))  # Fourier coefficients of η'
	        function dη( x :: Vector{Float64} )
	            return real.(exp.(1im*(x.-x₀)*k')*hatdη/length(k))
	        end

			# Define the Jacobian as a linear map
			# to be used with GMRES iterative solver
			function JacF( u )
				dx.= real.(ifft(u))
				fη .= fft(data.η(x+ϵ/ν*dx))
				δ  = meanf(fη)
				dF(φ) = φ+ϵ*sqrt(μ)*1im*Π⅔ .* ( cotanh(sqrt(μ)*(1+ϵ*δ)*k).*fft(dη(x+ϵ/ν*dx) .* ifft(φ))
								+ϵ/(1+ϵ*δ)* (xdcotanh(sqrt(μ)*(1+ϵ*δ)*k).*fη ) * mean(dη(x+ϵ/ν*dx) .* ifft(φ)) )
				return LinearMap(dF, length(u); issymmetric=false, ismutating=false)
			end

			# Define the Jacobian as a matrix
			# to be used with direct solver
			if method ==3
				FFT = exp.(-1im*k*(x.-x₀)');
				IFFT = exp.(1im*k*(x.-x₀)')/length(x);
				Id = Diagonal(ones(size(x)));
				Mean = ones(size(x))'/length(x)
				M(v) = Diagonal( v )

				function JacFMat( u )
					dx.=real.(ifft(u))
					fη .= fft(data.η(x+ϵ/ν*dx))
					δ  = meanf(fη)
		        	return Id + ϵ*sqrt(μ)*1im*Π⅔ .* (M(cotanh(sqrt(μ)*(1+ϵ*δ)*k)) * FFT * M(dη(x+ϵ/ν*dx)) * IFFT
								+ ϵ/(1+ϵ*δ)* (xdcotanh(sqrt(μ)*(1+ϵ*δ)*k).* fη )  * Mean * M(dη(x+ϵ/ν*dx)) * IFFT )
				end
			end

			# The iterative map to solve F(u)=u
			function iterate( u )
				if method==1  		# by contraction fix point algorithm
			    	return -F(u)
				elseif method==2	# by Newton algorithm with GMRES iterative solver to invert the Jacobian
					return gmres( JacF(u) , -F(u) ; reltol = tol/100, verbose=verbose )
				elseif method==3	# by Newton algorithm with direct solver to invert the Jacobian
					return - JacFMat( u ) \ F(u)
				else
					@error("In the function `WaterWaves`, the argument `method` should be 1, 2 or 3.")
				end

			end


			# initiate the iterative argument
			fη.=fft(data.η(x))
			δ = meanf(fη)
			u0 .= -ν*sqrt(μ)*1im*Π⅔ .*cotanh(sqrt(μ)*k*(1+ϵ*δ)).*fη

			# perform the iterative loop
			norm0=norm(fη)
			normdiff=norm0
			niter=0
			while  normdiff>tol*norm0 && niter<maxiter
				nu0 = iterate(u0)
				u0 += nu0
				normdiff=norm(nu0)
			    niter+=1
				if verbose
					@info "Relative error  $(normdiff/norm0)  at step  $niter"
				end
			end
			if verbose
				if niter == maxiter
				    @warn "The iterative solver did not converge. Maybe try a different method."
				    @warn "Estimated normalized error : $(normdiff/norm0)"
				else
					@info "The iterative solver converged in $niter iterations."
					@info "Estimated normalized error : $(normdiff/norm0)"
				end
			end

			# Constructs relevant variables from u0 the solution to F(u)=0
			z0 .= data.η(x+ϵ/ν*real.(ifft(u0)))
			v0 .= data.v(x+ϵ/ν*real.(ifft(u0))) .* (1 .+ ϵ/ν*real.(ifft(Π⅔.*∂ₓ.*u0)) )

			U=[ Π⅔ .* fft(z0)   Π⅔ .* fft(v0) ]  # Π⅔ for dealiasing
			U[ abs.(U).< ktol ].=0     # applies Krasny filter if ktol>0
 			return U
		  end
		end

		# Reconstruct physical variables from raw data
		# Return `(η,v)`, where
		# - `η` is the surface deformation;
		# - `v` is the derivative of the trace of the velocity potential.
		function mapfro(U)
		   	ξ  .= real.(sqrt(μ)*(1+ϵ*meanf(U[:,1]))*k)
	       	xv .= real.(-1im*sqrt(μ)*ifft( Π⅔ .*cotanh(ξ) .*  U[:,1] ))

		   return x + ϵ*xv, real.( ifft(U[:,1]) ) , real.( ifft(U[:,2]) )
		end

		# Water Waves equations are ∂t U = f(U)
		function f!(U)
		   	Dphi .= real.(ifft(U[:,2]))
			ξ .= sqrt(μ)*(1 .+ ϵ*meanf(U[:,1]))*k
			Dxv .= real.(sqrt(μ)*ifft( k .* cotanh(ξ) .* U[:,1]  ))
			Dz .= real.(ifft( ∂ₓ.*  U[:,1] ))

			J .= (1 .+ ϵ*Dxv).^2 + μ*(ϵ*Dz).^2
			M1 .= real.(-1im/sqrt(μ)*ifft( mytanh(ξ).* U[:,2] ))
			M2 .= real.( 1im*sqrt(μ)*ifft( cotanh(ξ) .* fft(M1./J) ))
			q0 = mean((1 .+ ϵ*Dxv).*M2 + ϵ*μ*Dz.*M1./J)

			U[:,2] .= -∂ₓ .* U[:,1] - ϵ* Π⅔ .* ∂ₓ .* fft( Dphi.*M2 + 1/2 .*(Dphi.^2 -μ*M1.^2)./J - q0*Dphi )/ν
			U[:,1] .=  Π⅔ .* fft( (1 .+ ϵ*Dxv).*M1./J - ϵ*Dz.*M2 + ϵ*q0*Dz )/ν

		end

		# Water waves equations are ∂t (U1,U2) = (f1(U1,U2) , f2(U1,U2))
		function f1!(U1,U2)
		   	Dphi .= real.(ifft(U2))
			ξ .= sqrt(μ)*(1 .+ ϵ*meanf(U1))*k
			Dxv .= real.(sqrt(μ)*ifft( k .* cotanh(ξ) .* U1 ))
			Dz .= real.(ifft( ∂ₓ.* U1))

			J .= (1 .+ ϵ*Dxv).^2 + μ*(ϵ*Dz).^2
			M1 = real.(-1im/sqrt(μ)*ifft( mytanh(ξ).* U2 ))
			M2 = real.( 1im*sqrt(μ)*ifft( cotanh(ξ) .* fft(M1./J) ))
			q0 = mean((1 .+ ϵ*Dxv).*M2 + ϵ*μ*Dz.*M1./J)

			U1 .= Π⅔ .* fft( (1 .+ ϵ*Dxv).*M1./J - ϵ*Dz.*M2 + ϵ*q0*Dz )/ν
		end
		function f2!(U1,U2)
		   	Dphi .= real.(ifft(U2))
			ξ .= sqrt(μ)*(1 .+ ϵ*meanf(U1))*k
			Dxv .= real.(sqrt(μ)*ifft( k .* cotanh(ξ) .* U1))
			Dz .= real.(ifft( ∂ₓ.* U1))

			J .= (1 .+ ϵ*Dxv).^2 + μ*(ϵ*Dz).^2
			M1 = real.(-1im/sqrt(μ)*ifft( mytanh(ξ).* U2 ))
			M2 = real.( 1im*sqrt(μ)*ifft( cotanh(ξ) .* fft(M1./J) ))
			q0 = mean((1 .+ ϵ*Dxv).*M2 + ϵ*μ*Dz.*M1./J)

			U2 .= -∂ₓ .* U1 - ϵ* Π⅔ .* ∂ₓ .* fft( Dphi.*M2 + 1/2 .*(Dphi.^2 -μ*M1.^2)./J - q0*Dphi )/ν
		end


		new(label, f!, f1!, f2!, mapto, mapfro, info )
    end
end
