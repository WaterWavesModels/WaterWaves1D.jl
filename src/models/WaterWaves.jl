export WaterWaves
#using Statistics #only for function mean in WaterWaves.jl...
using FFTW,LinearAlgebra

function cotanh( x :: Vector{Float64} )
	y=1 ./ tanh.(x+(x.==0))
	y[x.==0].=0
	return y
end
function mean(x)
	sum(x)/length(x)
end

"""
    WaterWaves(params)
	The Water Waves system, through conformal mapping

"""
mutable struct WaterWaves <: AbstractModel

    label   :: String
	mapto	:: Function
	mapfro	:: Function
	f!		:: Function

    function WaterWaves(param::NamedTuple)

		label = "water waves"
		datasize = 2
		μ 	= param.μ
		ϵ 	= param.ϵ
		mesh = Mesh(param)

		x 	= mesh.x
		k 	= mesh.k

    	∂ₓ	=  1im * mesh.k            # Differentiation
        Π⅔ 	= abs.(mesh.k) .< mesh.kmax * 2/3 # Dealiasing low-pass filter

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

		# Constructs conformal variables in the flat strip
		function mapto(data::InitialData)

			if norm(data.v(x))!=0
				error("non-zero initial data are not implemented yet. Sorry!")
			else

				function iter( u :: Vector{Float64} )
				    return imag.(ifft(sqrt(μ)*Π⅔ .*cotanh(sqrt(μ)*(1+ϵ*mean(data.η(x+ϵ*u)))*k).*fft(data.η(x+ϵ*u))))
				end

				# on démarre avec deux coefficients spécifiques:
				δ = mean(data.η(x))
				u0 = imag.(ifft(sqrt(μ)*Π⅔ .*cotanh(sqrt(μ)*k).*fft(data.η(x).-δ)))

				# détermination du point fixe de la fonction contractante:
				norm0=norm(data.η(x))
				tol_cont=1e-16
				max_iter=10000
				niter=0
				nu0=u0 .+ 1
				while norm(nu0-u0)>tol_cont*norm0 && niter<max_iter
				    nu0=u0
				    u0 = iter(nu0)
				    niter+=1
				end

				if niter==max_iter
				    @warn "The fix point algorithm did not converge"
				    @warn string("Estimated error : ",norm(nu0-u0)/norm0)
				else
					@warn string("The fix point algorithm converged in ",niter," iterations")
				end
				δ = mean(data.η(x+ϵ*u0))

				# Obtention des conditions initiales dans le domaine redressé:
				#x0 = x + ϵ*real.(ifft(Π⅔.*fft(u0))) #axe longitudinal
				z0 = δ .+ 1/sqrt(μ)*real.(ifft(1im*Π⅔.*tanh.(sqrt(μ)*k*(1+ϵ*δ)).*fft(u0)))
				phi0 = 0 .* z0

			end
			[z0 phi0] ######
		end

		# Reconstructs physical variables from conformal variables
		function mapfro(U)

				   ξ .= real.(sqrt(μ)*(1+ϵ*mean(U[:,1]))*k)
			       xv .= imag.(sqrt(μ)*ifft( cotanh(ξ) .* fft( U[:,1] )))

				   x + ϵ*xv, real.(U[:,1]) , real.(ifft(∂ₓ .* fft( U[:,2] )))
		end

		# Water Waves equations are ∂t U = f!(U)
		function f!(U)
			   	z .= U[:,1]
			   	phi .= U[:,2]
				ξ .= sqrt(μ)*(1 .+ ϵ*mean(z))*k
				xv .=imag.(sqrt(μ)*ifft( Π⅔ .* cotanh(ξ) .* fft( z )))
				Dxv .= real.(ifft(Π⅔ .* ∂ₓ.* fft(xv)))
				Dz .= real.(ifft(Π⅔ .* ∂ₓ.*fft(z)))
				Dphi .= real.(ifft(Π⅔ .* ∂ₓ.*fft(phi)))

				J .= (1 .+ ϵ*Dxv).^2 + μ*(ϵ*Dz).^2 # Jacobien
				M1 = imag.(1/sqrt(μ)*ifft(Π⅔ .* tanh.(ξ).* ∂ₓ .* fft(phi)))
				M2 = imag.( -sqrt(μ)*ifft(Π⅔ .* cotanh(ξ) .* fft(M1./J )))
				q0 = mean((1 .+ ϵ*Dxv).*M2 + ϵ*μ*Dz.*M1./J)

				U[:,1] .= (1 .+ ϵ*Dxv).*M1./J - ϵ.*Dz.*M2 + ϵ*q0*Dz
				U[:,2] .= -z -ϵ*Dphi.*M2 + 0.5*ϵ*μ*(M1.^2)./J - 0.5*ϵ*(Dphi.^2)./J + ϵ*q0*Dphi

			end


        new(label, mapto, mapfro, f!)
    end
end
