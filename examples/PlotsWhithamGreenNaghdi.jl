# #
# Reproduces the figures in the work of C. Klein and V. Duchêne
# on the Serre-Green-Naghdi and Whitham-Green-Naghdi systems
# #
#export PlotSolitaryWaveWGN1,PlotSolitaryWaveWGN2,PlotSolitaryWaveWGN3,PlotJacobianWGN,IntegrateSolitaryWaveWGN,StabilitySolitaryWaveWGN,IntegrateWGN

using ShallowWaterModels,FFTW,Plots,LinearAlgebra,ProgressMeter;gr()
#include("../src/dependencies.jl")

#---- Figures 1 and 2
"""
	`PlotSolitaryWaveWGN1(;kwargs)`

First method to compute the WGN solitary wave.
Uses GMRES-iterative method for the inversion of the Jacobian.

Arguments are all optional:
- `c` the velocity,
- `N` the number of collocation points,
- `L` the half-length of the mesh,
- `sav` a string used to save the figure as a .pdf.

Use with `PlotSolitaryWaveWGN1()` or e.g. `PlotSolitaryWaveWGN1(c=1.1,N=2^8,L=10*π)`
"""
function PlotSolitaryWaveWGN1(;c=2,N=2^10,L=10*π,sav=[])
	param = ( μ  = 1,
			ϵ  = 1,
        	N  = N,
            L  = L,
			c  = c
				)

	(η,u,v,mesh) = SolitaryWaveWhithamGreenNaghdi(
				param; SGN = false,
				method=2, α = 0,
				tol =  1e-16, max_iter=10,
				ktol =0*1e-11, gtol = 1e-16,
				iterative = true, q=1,
				verbose = false)

	(ηGN,uGN) = SolitaryWaveWhithamGreenNaghdi(
				param; SGN = true, max_iter=0)

	plt = plot(layout=(1,2))
	plot!(plt[1,1], mesh.x, [u uGN];
	  title=string("c=",c),
	  xlabel = "x",
	  ylabel = "u",
	  label=["WGN" "SGN"])
	plot!(plt[1,2], fftshift(mesh.k),
	  [log10.(abs.(fftshift(fft(u)))) log10.(abs.(fftshift(fft(uGN))))];
	  title="frequency",
	  label=["WGN" "SGN"])
	display(plt)
	if sav != [] savefig(string(sav,".pdf")); end
end


#---- Figure 3
"""
	PlotSolitaryWaveWGN2(;kwargs)

Second method to compute the WGN solitary wave. Same as `PlotSolitaryWaveWGN1`, but
uses non-iterative method for the inversion of the Jacobian.
To be used with higher values of the velocity (default is `c=20`).
"""
function PlotSolitaryWaveWGN2(;c=20,L=10*π,N=2^10,sav=[])
	param = ( μ  = 1,
		ϵ  = 1,
      	N  = N,
        L  = L,
		c = c )

	(η,u,v,mesh) = SolitaryWaveWhithamGreenNaghdi(
				param; SGN = false,
				method=2, α = 1, #ici α = 1 évite des oscillations important si c = 3 ou c = 20
				tol =  1e-14, max_iter=15,
				iterative = false, q=1,
				verbose = true)

	(ηGN,uGN) = SolitaryWaveWhithamGreenNaghdi(
				param; SGN = true, max_iter=0)

	plt = plot(layout=(1,2))
	plot!(plt[1,1], mesh.x, [u uGN];
	  title=string("c=",c),
	  xlabel = "u",
	  ylabel = "u/c",
	  label=["WGN" "SGN"])
	plot!(plt[1,2], fftshift(mesh.k),
	  [log10.(abs.(fftshift(fft(u)))) log10.(abs.(fftshift(fft(uGN))))];
	  title="frequency",
	  label=["WGN" "SGN"])
    display(plt)
	if sav != [] savefig(string(sav,".pdf")); end
end
#---- Figure 4
"""
	PlotSolitaryWaveWGN3(;kwargs)

Third method to compute the WGN solitary wave. Same as `PlotSolitaryWaveWGN1`, but
- uses non-iterative method for the inversion of the Jacobian; and
- uses a rescaled equation.
To be used with highest values of the velocity (default is `c=100`).
"""
function PlotSolitaryWaveWGN3(;c=100,L=10*π,N=2^10,sav=[])
	param = ( μ  = 1,
		ϵ  = 1,
      	N  = N,
        L  = L,
		c = c )

	(η,u,v,mesh) = SolitaryWaveWhithamGreenNaghdi(
				param; SGN = false,
				method=3, α = 1,
				tol =  1e-10, max_iter=10,
				iterative = false, q=1,
				verbose = true)
	(ηGN,uGN) = SolitaryWaveWhithamGreenNaghdi(
				param; SGN = true, max_iter=0)

	plt = plot(layout=(1,2))
	plot!(plt[1,1], mesh.x, [u/c uGN/c];
	  title=string("c=",c),
	  xlabel = "x",
	  ylabel = "u/c",
	  label=["WGN" "SGN"])
	plot!(plt[1,2], fftshift(mesh.k),
	  [log10.(abs.(fftshift(fft(u/c)))) log10.(abs.(fftshift(fft(uGN/c))))];
	  title="frequency",
	  label="WGN")
	display(plt)
	if sav != [] savefig(string(sav,".pdf")); end
end

#------ Figure 6
"""
	`PlotJacobianWGN(;kwargs)`

Computes the Jacobian matrix to be inverted in SGN or WGN equation.

Arguments are all optional:
- `c` the velocity,
- `N` the number of collocation points,
- `L` the half-length of the mesh,
- `SGN` uses SGN if `true`, and WGN if `false` (default is `false`),
- `sav` a string used to save the figure as a .pdf.

Use with `PlotJacobianWGN()` or e.g. `PlotJacobianWGN(c=20,N=2^10,L=10*π,SGN=true)`
"""
function PlotJacobianWGN(;c=20,L=10*π,N=2^10,SGN=false,sav=[])
	ϵ,μ,α=1,1,0

	(η,u,v,mesh) = SolitaryWaveWhithamGreenNaghdi(
				(c=c,L=L,N=N,μ=μ,ϵ=ϵ); SGN = SGN,
				method=2, α = 1, #ici α = 1 évite des oscillations important si c = 3 ou c = 20
				tol =  1e-12, max_iter=15,
				iterative = false, q=1,
				verbose = true)

	k,x=mesh.k,mesh.x
	if SGN == true
		F₀ = sqrt(μ)*1im * k
	else
		F₁ = tanh.(sqrt(μ)*abs.(k))./(sqrt(μ)*abs.(k))
		F₁[1] = 1
		F₀ = 1im * sqrt.(3*(1 ./F₁ .- 1)).*sign.(k)
	end
	x₀ = mesh.x[1]
	FFT = exp.(-1im*k*(x.-x₀)');
	IFFT = exp.(1im*k*(x.-x₀)')/length(x);
	M₀ = IFFT * Diagonal( F₀ )* FFT
	M(v) = Diagonal( v )


	dxu = real.(ifft(1im*k.*fft(u)))
	dxu ./= norm(dxu,2)
	hu = c ./(c .- ϵ*u)
	Fu = hu.* real.(ifft(F₀.*fft(u)))
	F2u = real.(ifft(F₀.*fft(hu.^2 .* Fu ) ))
	Du = c ./hu .+ 2*ϵ/3 * F2u ./ hu .+ ϵ^2/c * hu .* Fu.^2 .- (hu.^2)/c
	Jac = (-1/3 *M(1 ./ hu.^2 )* M₀ * M(hu.^3)* (c*M₀ .+ 3*ϵ * M( Fu ))
					.+ ϵ * M( hu .* Fu ) * M₀
					.+ M( Du ) .+ α*M( hu.^2 )*dxu*dxu' *M(1 ./ hu.^2 ) )
	Jacstar = -1/3 *M(1 ./ hu.^2 )* M₀ * M(hu.^3)* c*M₀

	plt = plot(layout=(1,2))
	surface!(plt[1,1],fftshift(k),fftshift(k)[N:-1:1],
		log10.(abs.(FFT*Jac*IFFT)),
		title = "Jacobian")
	surface!(plt[1,2],fftshift(k),fftshift(k)[N:-1:1],
		log10.(abs.(FFT*Jacstar*IFFT)),
		title = "non-diagonal part")
	display(plt)

	if sav != [] savefig(string(sav,".pdf")); end
end

#---- Figures 7 and 8
"""
	`IntegrateSolitaryWaveWGN(;kwargs)`

Integrates in time the SGN or WGN solitary wave.

Arguments are all optional:
- `c` the velocity,
- `N` the number of collocation points,
- `L` the half-length of the mesh,
- `T` final time of integration,
- `dt` timestep,
- `SGN` uses SGN if `true`, and WGN if `false` (default is `false`),
- `sav` a string used to save raw data and the figure as a .pdf.

Use with `IntegrateSolitaryWaveWGN()` or e.g. `IntegrateSolitaryWaveWGN(c=2,N=2^9,L=10*π,T=1,dt=1/2000,SGN=true)`
"""
function IntegrateSolitaryWaveWGN(;SGN=false,c=2,N=2^10,L=10*π,T=1,dt=1/2000,sav=[])
	if sav != [] ns=floor(Int,max(1,T/dt/100)) else ns=1 end

	param = ( μ  = 1, ϵ  = 1, c = c,
				N  = N,
	            L  = L,
	            T  = T,
	            dt = dt, ns=ns)
	if SGN == true
		(ηinit,uinit,vinit,mesh) = SolitaryWaveWhithamGreenNaghdi(
				param; SGN = true, max_iter=0)
		(ηfin,ufin,vfin) = SolitaryWaveWhithamGreenNaghdi(
				param; x₀ = c*T, SGN = true, max_iter=0)
	else
		(ηinit,uinit,vinit,mesh) = SolitaryWaveWhithamGreenNaghdi(
				param; method = 3, SGN = false, max_iter=10,α=1,verbose=true)
	  	(ηfin,ufin,vfin) = SolitaryWaveWhithamGreenNaghdi(
				param; x₀ = c*T, method = 3, SGN = false, max_iter=10,α=1,verbose=true)
	end
	init     = Init(mesh,ηinit,vinit)
	model = WhithamGreenNaghdi(param;SGN=SGN, ktol=0, iterate=true, precond = false)
	problem = Problem(model, init, param)
	solve!( problem )

	(ηcomp,vcomp,ucomp)   =  mapfrofull(model,last(problem.data.U))

	E(η,u,v) = sum(η.^2 .+ (1 .+ param.ϵ*η).*u.*v)
	# E1(η,u) = sum(η.^2 .+ (1 .+ param.ϵ*η).*u.^2 .+param.μ/3*(1 .+ param.ϵ*η).^3 .*(Dx(u).^2))
	# E2(η,u) = sum(η.^2 .+ (1 .+ param.ϵ*η).*u.^2 .-param.μ/3*u.*Dx((1 .+ param.ϵ*η).^3 .*(Dx(u))))
	dE(η1,u1,v1,η2,u2,v2) = sum(η1.^2-η2.^2) + sum((1 .+ param.ϵ*η1).*u1.*v1 - (1 .+ param.ϵ*η2).*u2.*v2)
	# dE1(η1,u1,η2,u2) = sum(η1.^2-η2.^2) + sum((1 .+ param.ϵ*η1).*u1.^2 - (1 .+ param.ϵ*η2).*u2.^2) -
	# 		param.μ/3*sum(u1.*Dx((1 .+ param.ϵ*η1).^3 .*(Dx(u1)))-u2.*Dx((1 .+ param.ϵ*η2).^3 .*(Dx(u2))))
	# dE2(η1,u1,η2,u2) = sum(η1.^2-η2.^2) + sum((1 .+ param.ϵ*η1).*u1.^2 - (1 .+ param.ϵ*η2).*u2.^2) +
	# 		param.μ/3*sum((1 .+ param.ϵ*η1).^3 .*(Dx(u1).^2) - (1 .+ param.ϵ*η2).^3 .*(Dx(u2).^2))

	print(string("normalized energy difference: ",dE(ηinit,uinit,vinit,ηcomp,ucomp,vcomp)/E(ηinit,uinit,vinit),"\n"))
	# print(string("normalized error: ",dE1(ηGN,uGN,ηfin,ufin)/E1(ηGN,uGN),"\n"))
	# print(string("normalized error: ",dE2(ηGN,uGN,ηfin,ufin)/E2(ηGN,uGN),"\n"))

	plt = plot(layout=(1,2))
	plot!(plt[1,1],fftshift(mesh.k),fftshift(log10.(abs.(fft(vinit))));
			title = "frequency",
			label = "initial",
			xlabel="x",
			ylabel="v")
	plot!(plt[1,1],fftshift(mesh.k),fftshift(log10.(abs.(fft(vcomp))));
			label = "final")
	plot!(plt[1,2],mesh.x,vfin-vcomp;
			title = string("error at time t=",problem.times.tfin),
			label="difference",
			xlabel="x",
			ylabel="v")
	display(plt)

	if sav != []
		savefig(string(sav,".pdf"));
		create_animation(problem;str=string(sav,"-anim.pdf"))
		plot_solution(problem)
		savefig(string(sav,"-final.pdf"));
		save(problem,sav);
	end
end

#------ Figures 9 to 11

"""
	`StabilitySolitaryWaveWGN(;kwargs)`

Integrates in time a pertubed SGN or WGN solitary wave.

Arguments are all optional:
- `p` (1,2,3 or 4) is the type of perturbation,
- `c` the velocity,
- `N` the number of collocation points,
- `L` the half-length of the mesh,
- `T` final time of integration,
- `dt` timestep,
- `SGN` uses SGN if `true`, and WGN if `false` (default is `false`),
- `sav` a string used to save raw data and the figure as a .pdf.

Use with `StabilitySolitaryWaveWGN()` or e.g. `IntegrateSolitaryWaveWGN(c=4,p=4,SGN=true)`

"""
function StabilitySolitaryWaveWGN(;p=2,c=2,N=2^10,L=10*π,T=10,dt=10/10^4,SGN=false,sav=[])
	if sav != [] ns=floor(Int,max(1,T/dt/100)) else ns=1 end
	if p == 1
		λ = 0.99
	elseif p == 2
		λ = 1.01
	elseif p == 3
		λ = -0.01
	elseif p == 4
		λ = +0.01
	else error("p must be in {1,2,3,4} ")
	end
	param = ( μ = 1, ϵ = 1, c = c, λ = λ,
				N  = N,
	            L  = L,
	            T  = T,
	            dt = dt, ns=ns,
				dealias = 0,
				SGN=SGN,
				ktol=0*1e-6,
				gtol=1e-12,
				iterate=true,
				precond = false)

	if SGN == true
		(η,u,v,mesh) = SolitaryWaveWhithamGreenNaghdi(
				param; SGN = true, max_iter=0)
	else
		(η,u,v,mesh) = SolitaryWaveWhithamGreenNaghdi(
				param; method = 3, SGN = false, max_iter=10,α=1,verbose=true)
	end

	if p == 1 || p == 2
		u= λ*u
	elseif p == 3 || p == 4
		u .+= λ*exp.(-mesh.x.^2)
	end
	k = mesh.k
	if SGN == true
		F₀ = sqrt(param.μ)*1im * k
	else
		F₁ = tanh.(sqrt(param.μ)*abs.(k))./(sqrt(param.μ)*abs.(k))
		F₁[1] = 1
		F₀ = 1im * sqrt.(3*(1 ./F₁ .- 1)).*sign.(k)
	end
	DxF(v) = real.(ifft(F₀ .* fft(v)))
	h = 1 .+ param.ϵ*η
	v = u - 1/3 ./h .* (DxF(h.^3 .*DxF(u)))

	init = Init(mesh,η,v)
	model = WhithamGreenNaghdi(param;SGN=param.SGN, dealias = param.dealias, ktol=param.ktol, gtol=param.gtol, iterate=param.iterate, precond = param.precond)
	problem = Problem(model, init, param)
	solve!( problem )

	if sav != [] save(problem,sav); end

	(ηfin,vfin,ufin)   =  mapfrofull(model,last(problem.data.U))


	# E(η) = sum(η.^2 ./ (1 .+ η) .+ param.μ/3*(1 .+ param.ϵ*η).^3 .*(Dx(η ./ (1 .+ η)).^2))
	E(η,u,v) = sum(η.^2 .+ (1 .+ param.ϵ*η).*u.*v)
	# E(η,u) = sum(η.^2 .+ (1 .+ param.ϵ*η).*u.^2 .+param.μ/3*(1 .+ param.ϵ*η).^3 .*(Dx(u).^2))
	dE(η1,u1,v1,η2,u2,v2) = sum(η1.^2-η2.^2) + sum((1 .+ param.ϵ*η1).*u1.*v1 - (1 .+ param.ϵ*η2).*u2.*v2)
	# dE(η1,u1,η2,u2) = sum(η1.^2-η2.^2) + sum((1 .+ param.ϵ*η1).*u1.^2 - (1 .+ param.ϵ*η2).*u2.^2) -
	# 		param.μ/3*sum(u1.*Dx((1 .+ param.ϵ*η1).^3 .*(Dx(u1)))-u2.*Dx((1 .+ param.ϵ*η2).^3 .*(Dx(u2))))

	print(string("normalized error: ",dE(η,u,v,ηfin,ufin,vfin)/E(η,u,v),"\n"))
	# print(string("normalized error: ",dE(ηGN,uGN,ηfin,ufin)/E(ηGN,uGN),"\n"))
	# print(string("normalized error: ",dE(ηinit,uinit,vinit,ηfin,ufin,vfin)/E(ηinit,uinit,vinit),"\n"))
	# print(string("normalized error: ",dE(ηinit,uinit,ηfin,ufin)/E(ηinit,uinit),"\n"))

	plt = plot(layout=(1,2))
	plot!(plt[1,1],fftshift(mesh.k),fftshift(log10.(abs.(fft(η))));
			title = "frequency",
			label = "initial")
	plot!(plt[1,1],fftshift(mesh.k),fftshift(log10.(abs.(fft(ηfin))));
			label = "final")
	plot!(plt[1,2],mesh.x,[η ηfin];
			title = string("at time t=",problem.times.tfin),
			label=["zeta initial" "zeta final"])
	plot!(plt[1,2],mesh.x,[u ufin];
			label=["u initial" "u final"])
	display(plt)
	if sav != [] savefig(string(sav,"-final.pdf")); end

	ts = problem.times.ts
	x = mesh.x
	k = mesh.k
	us=zeros(param.N,length(ts));
	# h = similar(Complex.(x))
	# fftu= Complex.(h)
	# Precond = Diagonal( 1 .+ param.μ*k.^2 )
	#
	# function compute(fftη,fftv)
	# 	h .= 1 .+ ifft(fftη)
	# 	function LL(hatu)
	# 		hatu- param.μ/3 *fft( 1 ./h .* ifft( 1im*k .* fft( h.^3 .* ifft( 1im*k .* hatu ) ) ) )
	# 	end
	# 	fftu.= gmres( LinearMap(LL, length(h); issymmetric=false, ismutating=false) , fftv ;
	# 				Pl = Precond )
	# 	return real.(ifft(fftu))
	# end
	@showprogress 1 for i in 1:length(ts)
		#us[:,i].=compute(problem.data.U[i][:,1],problem.data.U[i][:,2])
		us[:,i].=real.(ifft(problem.data.U[i][:,1]))
	end


	plt = plot()
	scatter!(plt,ts,maximum(abs.(us),dims=1)',
			title="maximum of surface deformation",
			label="",
			xlabel="time t")
	display(plt)

	if sav != []
		savefig(string(sav,"-norm.pdf"));
		plt=plot()
		my_cg = cgrad([:blue,:green])
		surface!(plt,x,ts,us',view_angle=(20,30), color = my_cg)
		display(plt)
		savefig(string(sav,"-evol.pdf"));

		create_animation(problem;str=string(sav,"-anim.pdf"))
	end
end

#---- Figures 14 to end
"""
	`IntegrateWGN(scenario;kwargs)

Integrates in time SGN or WGN with an initial data depending on a given `scenario`
- If `scenario = 1`, produces a dispersive shock wave for a unidirectional wave constructed from the Saint-Venant model.
- If `scenario = 2`, produces a dispersive shock wave for a unidirectional wave constructed from the Camassa-Holm model.
- If `scenario = 3`, the initial data challenges the non-cavitation assumption.

Other arguments are optional:
- `δ` is either
	- the square root of the shallowness parameter (the nonlinear parameter `ϵ=1`) if `scenario=1` or `2`;
	- the minimal depth if `scenario=3`
- `N` the number of collocation points,
- `L` the half-length of the mesh,
- `T` final time of integration,
- `dt` timestep,
- `dealias`: dealiasing with Orlicz rule `1-dealias/(dealias+2)` (default is `0`, i.e. no dealiasing);
- `SGN` uses SGN if `true`, and WGN if `false` (default is `false`),
- `sav` a string used to save raw data and the figure as a .pdf.

Use with `IntegrateWGN(s)` or e.g. `IntegrateWGN(3;δ=0.1,SGN=true,N=2^12,L=3*π,T=0.6,dt=0.6/10^4,dealias=1)`

"""
function IntegrateWGN(scenario;δ=0.1,N=2^11,L=3*π,T= 5,dt = 5/10^4,SGN=false,dealias=0,sav=[])

	if sav != [] ns=floor(Int,max(1,T/dt/100)) else ns=1 end
	if scenario == 1 || scenario == 2
		μ  = δ^2
	else
		μ = 1
	end
	param = ( μ  = μ, ϵ  = 1,
				N  = N, L  = L,
	            T  = T, dt = dt,
				ns=ns )

	mesh=Mesh(param)
	if scenario == 1
		η= exp.(-(mesh.x .+3).^2)
		u= 2*sqrt.(1 .+ param.ϵ*η) .-2
	elseif scenario == 2
		function krasny!(v)
			v[abs.(v).<1e-14].=0
			return v
		end
		Dx(v) = real.(ifft( 1im*mesh.k.* krasny!(fft(v))))
		Dx2(v) = Dx(Dx(v))
		ϵ = param.ϵ;μ=param.μ;
		w = - (mesh.x).* exp.(-(mesh.x).^2)
		u = w .+ μ/12 *Dx2(w) .+ μ*ϵ/6* w.*Dx2(w)
		η = u .+ ϵ/4*u.^2 .- μ/6* Dx2( u+3*ϵ/4* u.^2) .- μ*ϵ/6 * u .* Dx2(u) .- 5*μ*ϵ/48 * Dx(u).^2
	elseif scenario == 3
		η = -(1-δ)*exp.(-mesh.x.^2)
		u = -mesh.x.*exp.(-mesh.x.^2)
	else
		error("the first argument must be 1, 2 or 3")
	end
	if SGN == true
		F₀ = sqrt(param.μ)*1im*mesh.k
	else
		F₁ 	= tanh.(sqrt(param.μ)*abs.(mesh.k))./(sqrt(param.μ)*abs.(mesh.k))
		F₁[1]= 1                 # Differentiation
		F₀   = 1im * sqrt.(3*(1 ./F₁ .- 1)).*sign.(mesh.k)
	end
	DxF(v) = real.(ifft(F₀ .* fft(v)))
	h=1 .+ param.ϵ*η
	v= u - 1/3 ./h .* (DxF(h.^3 .*DxF(u)))

	init     = Init(mesh,η,v)
	model = WhithamGreenNaghdi(param;SGN=SGN, ktol=0, gtol=1e-12, iterate=true, precond = false, dealias = dealias)
	problem = Problem(model, init, param)
	solve!( problem )

	(ηfin,vfin,ufin)   =  mapfrofull(model,last(problem.data.U))
	E(η,u,v) = sum(η.^2 .+ (1 .+ param.ϵ*η).*u.*v)
	dE(η1,u1,v1,η2,u2,v2) = sum(η1.^2-η2.^2) + sum((1 .+ param.ϵ*η1).*u1.*v1 - (1 .+ param.ϵ*η2).*u2.*v2)

	print(string("normalized error: ",dE(η,u,v,ηfin,ufin,vfin)/E(η,u,v),"\n"))

	fftηfin=last(problem.data.U)[:,1]

	plt = plot(layout=(1,2))
	plot!(plt[1,1],fftshift(mesh.k),fftshift(log10.(abs.(fftηfin)));
			title = "frequency",label="")
	plot!(plt[1,2],mesh.x,real.(ifft(fftηfin));
			title = string("surface deformation at time t=",problem.times.tfin),
			label="")
	if scenario == 1
		xlims!(plt[1,2],(0,8))
	end
	display(plt)

	if sav != []
		savefig(string(sav,".pdf"));
		p = plot(layout=(2,1))
		save(problem,sav);
		ts = problem.times.ts
		x = mesh.x
		k = mesh.k
		us=zeros(param.N,length(ts));
		@showprogress 1 for i in 1:length(ts)
			us[:,i].=real.(ifft(problem.data.U[i][:,1]))
		end
		plt=plot()
		my_cg = cgrad([:blue,:green])
		surface!(plt,x,ts,us',view_angle=(20,30), color = my_cg)
		display(plt)
		savefig(string(sav,"-evol.pdf"));
		anim = @animate for i in range(1,length(ts))
			plt = plot(layout=(1,2))
			fftη=problem.data.U[i][:,1]
			plot!(plt[1,1],fftshift(mesh.k),fftshift(log10.(abs.(fftη)));
					title = "frequency",
					label = "")
			plot!(plt[1,2],mesh.x,real.(ifft(fftη));
					title = string("surface at time t=",problem.times.ts[i]),
					label="")
		end

		gif(anim, string(sav,"-anim.gif"), fps=15); nothing
	end
end
test
