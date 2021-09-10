export WaterWaves
#using Statistics #only for function mean in WaterWaves.jl...
using FFTW,LinearAlgebra

"""
    WaterWaves(params)

Define an object of type `AbstractModel` in view of solving the water waves system
(via conformal mapping).

# Argument
`param` is of type `NamedTuple` and must contain
- dimensionless parameters `ϵ` (nonlinearity) and `μ` (dispersion);
- numerical parameters to construct the mesh of collocation points as `mesh = Mesh(param)`.

# Return values
Generate necessary ingredients for solving an initial-value problem via `solve!` and in particular
1. a function `WaterWaves.f!` to be called in the time-integration solver;
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
	param	:: NamedTuple
	kwargs  :: NamedTuple


    function WaterWaves(param::NamedTuple;dealias=0,method=1,tol=1e-16,maxiter=100,ktol=0,verbose=true)

		# Preparation
		label = "water waves"
		kwargs = (dealias=dealias,tol=tol,maxiter=maxiter,verbose=verbose)

		μ 	= param.μ
		ϵ 	= param.ϵ

		mesh = Mesh(param)
		x 	= mesh.x
		x₀ = x[1]
		k 	= mesh.k

    	∂ₓ	=  1im * mesh.k             # Differentiation Fourier multiplier
		K = mesh.kmax * (1-dealias/(2+dealias))
		Π⅔ 	= abs.(mesh.k) .<= K 		# Dealiasing low-pass filter
		if dealias == 0
			if verbose @info "no dealiasing" end
			Π⅔ 	= ones(size(mesh.k))
		elseif verbose
			@info string("dealiasing : spectral scheme for power ", dealias + 1," nonlinearity ")
		end

		# preallocate some variables for memory use
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
		q0 = 0

		# some useful functions
		function cotanh( x )
			y=1 ./ tanh.(x+(x.==0))
			y[x.==0].=0
			return y
		end
		function xdcotanh( x )
			return -x ./ (sinh.(x+(x.==0)).^2)
		end
		function meanf(f) #computes the mean of a function, being given its FFT
			real(f[1])/length(f)
		end
		function mean(f) #computes the mean of a function
			sum(f)/length(f)
		end


		# Construct conformal variables in the flat strip from initial data
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
				fη .= fft(data.η(x+ϵ*dx))
			    return u+sqrt(μ)*1im*Π⅔ .*cotanh(sqrt(μ)*(1+ϵ*meanf(fη))*k).*fη
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
				fη .= fft(data.η(x+ϵ*dx))
				δ  = meanf(fη)
				dF(φ) = φ+ϵ*sqrt(μ)*1im*Π⅔ .* ( cotanh(sqrt(μ)*(1+ϵ*δ)*k).*fft(dη(x+ϵ*dx) .* ifft(φ))
								+ϵ/(1+ϵ*δ)* (xdcotanh(sqrt(μ)*(1+ϵ*δ)*k).*fη ) * mean(dη(x+ϵ*dx) .* ifft(φ)) )
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
					fη .= fft(data.η(x+ϵ*dx))
					δ  = meanf(fη)
		        	return Id + ϵ*sqrt(μ)*1im*Π⅔ .* (M(cotanh(sqrt(μ)*(1+ϵ*δ)*k)) * FFT * M(dη(x+ϵ*dx)) * IFFT
								+ ϵ/(1+ϵ*δ)* (xdcotanh(sqrt(μ)*(1+ϵ*δ)*k).* fη )  * Mean * M(dη(x+ϵ*dx)) * IFFT )
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
			u0 .= -sqrt(μ)*1im*Π⅔ .*cotanh(sqrt(μ)*k*(1+ϵ*δ)).*fη

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
					@info string("Relative error ", normdiff/norm0, " at step ", niter)
				end
			end
			if verbose
				if niter == maxiter
				    @warn "The iterative solver did not converge. Maybe try a different method."
				    @warn string("Estimated normalized error : ",normdiff/norm0)
				else
					@info string("The iterative solver converged in ",niter," iterations.")
					@info string("Estimated normalized error : ",normdiff/norm0)
				end
			end

			# Constructs relevant variables from u0 the solution to F(u)=0
			z0 .= data.η(x+ϵ*real.(ifft(u0)))
			v0 .= data.v(x+ϵ*real.(ifft(u0))) .* (1 .+ ϵ*real.(ifft(Π⅔.*∂ₓ.*u0)) )

			U=[ Π⅔ .* fft(z0)   Π⅔ .* fft(v0) ]  # Π⅔ for dealiasing
			U[ abs.(U).< ktol ].=0     # applies Krasny filter if ktol>0
 			return U
		  end
		end

		# Reconstructs physical variables from conformal variables
		function mapfro(U)
			z  .= real.(ifft(U[:,1]))
		   	ξ  .= real.(sqrt(μ)*(1+ϵ*mean(z))*k)
	       	xv .= real.(-1im*sqrt(μ)*ifft( Π⅔ .*cotanh(ξ) .*  U[:,1] ))

		   return x + ϵ*xv, z , real.( ifft(U[:,2]) )
		end

		# Water Waves equations are ∂t U = f!(U)
		function f!(U)
		   	Dphi .= real.(ifft(U[:,2]))
			ξ .= sqrt(μ)*(1 .+ ϵ*meanf(U[:,1]))*k
			Dxv .= real.(sqrt(μ)*ifft( k .* cotanh(ξ) .* U[:,1]  ))
			Dz .= real.(ifft( ∂ₓ.*  U[:,1] ))

			J .= (1 .+ ϵ*Dxv).^2 + μ*(ϵ*Dz).^2
			M1 = real.(-1im/sqrt(μ)*ifft( tanh.(ξ).* U[:,2] ))
			M2 = real.( 1im*sqrt(μ)*ifft( cotanh(ξ) .* fft(M1./J) ))
			q0 = mean((1 .+ ϵ*Dxv).*M2 + ϵ*μ*Dz.*M1./J)

			U[:,2] .= -∂ₓ .* U[:,1] - ϵ* Π⅔ .* ∂ₓ .* fft( Dphi.*M2 + 1/2 .*(Dphi.^2 -μ*M1.^2)./J - q0*Dphi )
			U[:,1] .=  Π⅔ .* fft( (1 .+ ϵ*Dxv).*M1./J - ϵ*Dz.*M2 + ϵ*q0*Dz )

		end

		# Define (f1!,f2!) a splitting of the function f! used in symplectic schemes
		function f1!(U1,U2)
		   	Dphi .= real.(ifft(U2))
			ξ .= sqrt(μ)*(1 .+ ϵ*meanf(U1))*k
			Dxv .= real.(sqrt(μ)*ifft( k .* cotanh(ξ) .* U1 ))
			Dz .= real.(ifft( ∂ₓ.* U1))

			J .= (1 .+ ϵ*Dxv).^2 + μ*(ϵ*Dz).^2
			M1 = real.(-1im/sqrt(μ)*ifft( tanh.(ξ).* U2 ))
			M2 = real.( 1im*sqrt(μ)*ifft( cotanh(ξ) .* fft(M1./J) ))
			q0 = mean((1 .+ ϵ*Dxv).*M2 + ϵ*μ*Dz.*M1./J)

			U1 .= Π⅔ .* fft( (1 .+ ϵ*Dxv).*M1./J - ϵ*Dz.*M2 + ϵ*q0*Dz )
		end
		function f2!(U1,U2)
		   	Dphi .= real.(ifft(U2))
			ξ .= sqrt(μ)*(1 .+ ϵ*meanf(U1))*k
			Dxv .= real.(sqrt(μ)*ifft( k .* cotanh(ξ) .* U1))
			Dz .= real.(ifft( ∂ₓ.* U1))

			J .= (1 .+ ϵ*Dxv).^2 + μ*(ϵ*Dz).^2
			M1 = real.(-1im/sqrt(μ)*ifft( tanh.(ξ).* U2 ))
			M2 = real.( 1im*sqrt(μ)*ifft( cotanh(ξ) .* fft(M1./J) ))
			q0 = mean((1 .+ ϵ*Dxv).*M2 + ϵ*μ*Dz.*M1./J)

			U2 .= -∂ₓ .* 1 - ϵ* Π⅔ .* ∂ₓ .* fft( Dphi.*M2 + 1/2 .*(Dphi.^2 -μ*M1.^2)./J - q0*Dphi )
		end


		new(label, f!, f1!, f2!, mapto, mapfro, param, kwargs )
    end
end
