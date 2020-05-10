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
	x   	:: Vector{Float64}
	k   	:: Vector{Float64}
    ∂ₓ      :: Vector{Complex{Float64}}
    Π⅔      :: BitArray{1}
	μ 		:: Float64
    ϵ 		:: Float64
	z    	:: Vector{Float64}
	phi    	:: Vector{Float64}
	ξ    	:: Vector{Float64}
	xv    	:: Vector{Float64}
	Dxv   	:: Vector{Float64}
	Dz    	:: Vector{Float64}
	Dphi    :: Vector{Float64}
	J    	:: Vector{Float64}
	M1    	:: Vector{Float64}
	M2    	:: Vector{Float64}
	q0 		:: Float64

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

		"""
		    mapto(WaterWaves, data)
			Constructs conformal variables in the flat strip

		"""
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

		"""
		    mapfro(WaterWaves, data)
			Reconstructs physical variables from conformal variables

		"""
		function mapfro(U)

				   ξ .= real.(sqrt(μ)*(1+ϵ*mean(U[:,1]))*k)
			       xv .= imag.(sqrt(μ)*ifft( cotanh(ξ) .* fft( U[:,1] )))

				   x + ϵ*xv, real.(U[:,1]) , real.(ifft(∂ₓ .* fft( U[:,2] )))
		end

		function f!(U)			# utile uniquement pour solve2!
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


        new(label, mapto, mapfro, f!, x, k, ∂ₓ, Π⅔, μ, ϵ, z, phi, ξ, xv, Dxv, Dz, Dphi, J, M1, M2, q0)
    end
end

function (m::WaterWaves)(U)	# utile uniquement pour solve!


   	m.z .= U[:,1]
   	m.phi .= U[:,2]
	m.ξ .= sqrt(m.μ)*(1 .+ m.ϵ*mean(m.z))*m.k
	m.xv .=imag.(sqrt(m.μ)*ifft( m.Π⅔ .* cotanh(m.ξ) .* fft( m.z )))
	m.Dxv .= real.(ifft(m.Π⅔ .* m.∂ₓ.* fft(m.xv)))
	m.Dz .= real.(ifft(m.Π⅔ .* m.∂ₓ.*fft(m.z)))
	m.Dphi .= real.(ifft(m.Π⅔ .* m.∂ₓ.*fft(m.phi)))

	m.J .= (1 .+ m.ϵ*m.Dxv).^2 + m.μ*(m.ϵ*m.Dz).^2 # Jacobien
	m.M1 = imag.(1/sqrt(m.μ)*ifft(m.Π⅔ .* tanh.(m.ξ).* m.∂ₓ .* fft(m.phi)))
	m.M2 = imag.( -sqrt(m.μ)*ifft(m.Π⅔ .* cotanh(m.ξ) .* fft(m.M1./m.J )))
	m.q0 = mean((1 .+ m.ϵ*m.Dxv).*m.M2 + m.ϵ*m.μ*m.Dz.*m.M1./m.J)

	U[:,1] .= (1 .+ m.ϵ*m.Dxv).*m.M1./m.J - m.ϵ.*m.Dz.*m.M2 + m.ϵ*m.q0*m.Dz
	U[:,2] .= -m.z -m.ϵ*m.Dphi.*m.M2 + 0.5*m.ϵ*m.μ*(m.M1.^2)./m.J - 0.5*m.ϵ*(m.Dphi.^2)./m.J + m.ϵ*m.q0*m.Dphi


end
