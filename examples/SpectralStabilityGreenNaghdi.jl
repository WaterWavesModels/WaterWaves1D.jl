# #
# Study of the spectral stability of the Serre-Green-Naghdi equations.
# #
@info "Define functions Spectrum,figspecCW,figspecSW,Experiments"

using WaterWaves1D,LinearAlgebra,FFTW,ProgressMeter,Plots;
include("../src/initialdata/SolitaryWaveWhithamGreenNaghdi.jl")
include("../src/initialdata/SolitaryWaveSerreGreenNaghdi.jl")
include("../src/initialdata/CnoidalWaveSerreGreenNaghdi.jl")

"""
    `Spectrum(param,η,u,c,ν;kwargs)`

Computes the spectrum about a specified travelling solution of the
the linearized Green-Naghdi model.

# Argument
- `param` should be a NamedTuple
	e.g. `param = ( μ  = 0.1, ϵ  = 1, N  = 2^9, L  = 10*π)`

- `(η,u)` are real vectors of the aforementioned solution where
where

    - `η` is the surface deformation;
    - `u` corresponds to the layer-averaged velocity.
- 'c' is the velocity of the solution,
- 'ν' is a Floquet coefficient.

## Keywords
- `SGN :: Bool`: if `true` computes the Serre-Green-Naghdi (SGN) instead of Whitham-Green-Naghdi (WGN) system (default is `false`);

# Return values
an Eigen factorization object `L` which
contains the eigenvalues in `L.values` and the
eigenvectors in the columns of the matrix `L.vectors`.
(The `k`th eigenvector can be obtained from the slice `F.vectors[:, k]``.)

"""
function Spectrum(param,η,u,c,ν)

	μ 	= param.μ
	ϵ 	= param.ϵ
	mesh = Mesh(param)
	k = mesh.k
	x 	= mesh.x
	x₀ = mesh.x[1]
	∂ₓ	=  1im * mesh.k
	F₁ = 1 ./(1 .+ μ/3*k.^2)
	F₀ = sqrt(μ)*∂ₓ
	FFT = exp.(-1im*k*(x.-x₀)');
	IFFT = exp.(1im*k*(x.-x₀)')/length(x);
	Id = Diagonal(ones(size(x)));


	h = 1 .+ ϵ*η
	DxF(v) = ifft(F₀ .* fft(v))
	v= u - 1/3 ./h .* (DxF(h.^3 .*DxF(u)))

	Dx = IFFT*Diagonal(∂ₓ)*FFT + ν*Id#Diagonal(m.x./sqrt.(1 .+ m.x.^2))
	DxF₀ = IFFT*Diagonal(F₀)*FFT + ν*Id#Diagonal(m.x./sqrt.(1 .+ m.x.^2))

	A11 = -Dx*Diagonal(ϵ*u)
	A12 = -Dx*Diagonal(h)
	A13 = zeros(size(A11))
	A21 = -Dx* ( Id - Diagonal( h.*(DxF(ϵ*u).^2) ) )
	A22 = -Dx*( ϵ*Diagonal( v - u ) - Diagonal( (h.^2).*DxF(ϵ*u)) *DxF₀ )
	A23 = -Dx*Diagonal(ϵ*u)

	B1 = -Diagonal(1 ./h)*DxF₀*Diagonal( (h.^2).*DxF(ϵ*u) ) + 1/3*Diagonal(DxF(h.^3 .* DxF(ϵ*u))./h.^2)
	B2 = Id - 1/3*Diagonal(1 ./h)*DxF₀*Diagonal(h.^3)*DxF₀

	C2 = inv(B2)
	C1 = -C2*B1

	M11 = A11+A12*C1 + c*Dx
	M12 = A13+A12*C2
	M21 = A21+A22*C1
	M22 = A23+A22*C2 + c*Dx

	return eigen([ M11 M12  ; M21 M22 ])
end

function ell(a₀,a₁,k)
	h₀ = a₀
	h₂ = a₁ + a₀
	h₁ = (1-k^2)*a₁ + a₀
	return (h₀=h₀,h₁=h₁,h₂=h₂)
end

"""
	`figspecCW(p;N,P,SGN)`

Spectrum of the linearized SGN equations about a cnoidal wave.

- Argument `p` is a NamedTuple and should contain `h₀<h₁<h₂` definig the cnoidal wave.
- `N` is the number of collocation points.
- `P` is the number of periods in the mesh (default is `P=1`).
- `M` is the number of (equally spaced) Floquet points conputed.
- Solves SGN if `SGN` is `true`, WGN otherwise (not implemented yet).

"""
function figspecCW(p;N=2^6,P=1,M=2^6,SGN=true::Bool)
	param = merge( p,( μ  = 1, ϵ  = 1, N  = N, P = P ))

	(η,u,v,mesh,parc)=CnoidalWaveSerreGreenNaghdi(param;P=P)
	@show parc

	plt = plot(layout=(2,1))
	plot!(plt[1,1], mesh.x, [η u],label=["η" "u"])
	plot!(plt[2,1], fftshift(mesh.k),
	  [abs.(fftshift(fft(η))).+eps() abs.(fftshift(fft(u))).+eps()],
	  label=["|η̂|" "|û|"],yscale=:log10)
	display(plt)

	param2=merge(param,(L=P*parc.λ,))

	σ = Spectrum(param2,η,u,parc.c,1im*0.).values
	val = σ#[Int(7*end/8):end]
	@showprogress 1 for ν=-π/parc.λ/P/2:π/parc.λ/P/M:π/parc.λ/P/2
		σ = Spectrum(param2,η,u,parc.c,1im*ν).values
		val=[val;σ]
		# uncomment the following lines to stop when a unstable mode has been found
		# if maximum(real.(σ))>0.1
		# 	return Spectrum(param,η,u,parc.c,1im*ν),merge(parc,(ν=ν,)),η,u,v,mesh
		# 	break
		# end
	end
	#plt = scatter(values,label="")
	#display(plt)
	return val,parc,η,u,v,mesh

end

"""
	figspecSW(c;N,L,SGN)

Spectrum of the linearized SGN equations about a cnoidal wave.

- Argument `c` is the velocity definig the solitary wave.
- `N` is the number of collocation points.
- `L` is half-length of the mesh.
- `M` is the number of (equally spaced) Floquet points conputed.
- Solves SGN if `SGN` is `true`, WGN otherwise.

"""
function figspecSW(c;N=2^6,L=10*π,M=2^6,SGN=true::Bool)
	param = ( c=c, μ  = 1, ϵ  = 1, N  = N, L = L )

	if SGN == true
		(η,u,v,mesh)=SolitaryWaveSerreGreenNaghdi(param)
	else
		(η,u,v,mesh)=SolitaryWaveWhithamGreenNaghdi(param;method=3,verbose=true)
	end
	plt = plot(layout=(2,1))
	plot!(plt[1,1], mesh.x, [η u],label=["η" "u"])
	plot!(plt[2,1], fftshift(mesh.k),
	  [abs.(fftshift(fft(η))).+eps() abs.(fftshift(fft(u))).+eps()],
	  label=["|η̂|" "|û|"],yscale=:log10)
	display(plt)


	σ = Spectrum(param,η,u,c,1im*0.).values
	val = σ#[Int(7*end/8):end]
	@showprogress 1 for ν=-π/L/2:π/L/M:π/L/2
		σ = Spectrum(param,η,u,c,1im*ν).values
		val=[val;σ]
		# uncomment the following lines to stop when a unstable mode has been found
		# if maximum(real.(σ))>0.1
		# 	return Spectrum(param,η,u,c,1im*ν),η,u,v,mesh
		# 	break
		# end
	end
	#plt = scatter(values,label="")
	#display(plt)
	return val,η,u,v,mesh

end

"""
	Experiments(exp)

Some numerical experiments, plotting the spectrum of the linearized SGN equations about various travelling waves
- `exp=1` : a cnoidal wave with very small amplitude
- `exp=2` : a cnoidal wave with very large period
- `exp=3` : several cnoidal waves with increasing periods
- `exp=4` : a solitary wave

Return the plot.
"""
function Experiments(exp)

if exp == 1 # experiment 1: very small amplitude
	σ,param,η,u,v,mesh=figspecCW(ell(0.3,0.005,0.75);P=1,N=2^5)
	plt=scatter(σ)
	ylims!(0,0.1)
	σ,param,η,u,v,mesh=figspecCW(ell(0.3,0.005,0.75);P=1,N=2^5,M=2^19) # on peut limiter les Floquet dans (-π/parc.λ/P/20,π/parc.λ/P/20)
	scatter(σ[abs.(real.(σ)) .> 1e-10])
	#ylims!(0.03158,0.0316)
	ylims!(0.0113725,0.011375)

elseif exp == 2 # experiment 2: very large period
	σ,param,η,u,v,mesh=figspecCW(ell(0.3,0.1,0.999999999999999);P=1,N=2^8,M=2^8)
	plt=scatter(σ)
	ylims!(0.2,0.25)
	ylims!(0.23,0.235)

elseif exp == 3 # experiment 3 : growing period
	v1,param,η,u,v,mesh=figspecCW((h₀=0.99,h₁=1,h₂=11);P=1,N=2^7)
	v2,param,η,u,v,mesh=figspecCW((h₀=0.999999,h₁=1,h₂=11);P=1,N=2^7)
	v3,param,η,u,v,mesh=figspecCW((h₀=0.9999999999,h₁=1,h₂=11);P=1,N=2^7)
	v4,param,η,u,v,mesh=figspecCW((h₀=0.99999999999999,h₁=1,h₂=11);P=1,N=2^7)

	plt=scatter(v1[abs.(v1).<15],label="a=10,L=6")
	scatter!(plt,v2[abs.(v2).<15],label="a=10,L=12")
	scatter!(plt,v3[abs.(v3).<15],label="a=10,L=18")
	scatter!(plt,v4[abs.(v4).<15],label="a=10,L=22")

elseif exp == 4# experiment 4 : solitary wave
	σ,η,u,v,mesh=figspecSW(3;L=10*π,N=2^8)
	plt=scatter(σ)
end
display(plt)
return plt
end
