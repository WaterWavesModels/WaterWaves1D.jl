export figure7,figure8,figure10,figure12

"""
Reproduces some figures in the work of C. Klein and V. Duchêne
"""
#using ShallowWaterModels
include("../src/dependencies.jl")


#---- Figure 7

"""
	`figure7(c;sav)`

SGN solitary wave integrated in time.
- Argument 'c' is the velocity of the wave
- Optional argument 'sav': if a string is given, then saves data and figures.

"""
function figure7(c::Real;sav=[])
	param = ( μ  = 1, ϵ  = 1, c=c,
				N  = 2^9,
	            L  = 10*π,
	            T  = 1,
	            dt = 1/2000, ns=1)
	function solη(x,c,ϵ,μ)
		(c^2-1)/ϵ*sech(sqrt(3*(c^2-1)/(c^2)/μ)/2*x)^2
	end
	function solu(x,c,ϵ,μ)
		c*solη(x,c,ϵ,μ) / (1 + ϵ * solη(x,c,ϵ,μ))
	end
	mesh=Mesh(param)
	ηGN= solη.(mesh.x,c,param.ϵ,param.μ)
	hGN=(1 .+ param.ϵ*ηGN)
	uGN= solu.(mesh.x,c,param.ϵ,param.μ)
	Dx(v) = real.(ifft(1im*mesh.k .* fft(v)))
	vGN= uGN - param.μ/3 ./hGN .* (Dx(hGN.^3 .*Dx(uGN)))

	init     = Init(mesh,ηGN,vGN)
	model = WhithamGreenNaghdi(param;SGN=true, ktol=0, iterate=true, precond = false)
	problem = Problem(model, init, param)

	@time solve!( problem )

	(ηfin,vfin,ufin)   =  mapfro(model,last(problem.data.U))
	(ηinit,vinit,uinit) = mapfro(model,first(problem.data.U))


	E(η,u,v) = sum(η.^2 .+ (1 .+ param.ϵ*η).*u.*v)
	E1(η,u) = sum(η.^2 .+ (1 .+ param.ϵ*η).*u.^2 .+param.μ/3*(1 .+ param.ϵ*η).^3 .*(Dx(u).^2))
	E2(η,u) = sum(η.^2 .+ (1 .+ param.ϵ*η).*u.^2 .-param.μ/3*u.*Dx((1 .+ param.ϵ*η).^3 .*(Dx(u))))
	dE(η1,u1,v1,η2,u2,v2) = sum(η1.^2-η2.^2) + sum((1 .+ param.ϵ*η1).*u1.*v1 - (1 .+ param.ϵ*η2).*u2.*v2)
	dE1(η1,u1,η2,u2) = sum(η1.^2-η2.^2) + sum((1 .+ param.ϵ*η1).*u1.^2 - (1 .+ param.ϵ*η2).*u2.^2) -
			param.μ/3*sum(u1.*Dx((1 .+ param.ϵ*η1).^3 .*(Dx(u1)))-u2.*Dx((1 .+ param.ϵ*η2).^3 .*(Dx(u2))))
	dE2(η1,u1,η2,u2) = sum(η1.^2-η2.^2) + sum((1 .+ param.ϵ*η1).*u1.^2 - (1 .+ param.ϵ*η2).*u2.^2) +
			param.μ/3*sum((1 .+ param.ϵ*η1).^3 .*(Dx(u1).^2) - (1 .+ param.ϵ*η2).^3 .*(Dx(u2).^2))

	print(string("normalized error: ",dE(ηGN,uGN,vGN,ηfin,ufin,vfin)/E(ηGN,uGN,vGN),"\n"))
	print(string("normalized error: ",dE1(ηGN,uGN,ηfin,ufin)/E1(ηGN,uGN),"\n"))
	print(string("normalized error: ",dE2(ηGN,uGN,ηfin,ufin)/E2(ηGN,uGN),"\n"))

	print(string("normalized error: ",dE(ηinit,uinit,vinit,ηfin,ufin,vfin)/E(ηinit,uinit,vinit),"\n"))
	print(string("normalized error: ",dE1(ηinit,uinit,ηfin,ufin)/E1(ηinit,uinit),"\n"))
	print(string("normalized error: ",dE2(ηinit,uinit,ηfin,ufin)/E2(ηinit,uinit),"\n"))



	plt = plot(layout=(1,2))
	plot!(plt[1,1],fftshift(mesh.k),fftshift(log10.(abs.(fft(vinit))));
			title = "frequency",
			label = "initial")
	plot!(plt[1,1],fftshift(mesh.k),fftshift(log10.(abs.(fft(vfin))));
			label = "final")
	plot!(plt[1,2],mesh.x,ufin-solu.(mesh.x.-c*param.T,c,param.ϵ,param.μ);
			title = string("error at time t=",problem.times.tfin),
			label="difference")
	display(plt)

	if sav != []
		savefig(string("fig7-",sav,".pdf"));
		#p = plot(layout=(2,1))
		#fig_problem!(p,problem)
		#savefig(string("fig7a-",sav,".pdf"));
		save(problem,string("fig7-",sav,".pdf"));
	end
end
#figure7(c) # c is the velocity
#------ Figure 8

"""
	figure8(c;sav)

WGN solitary wave integrated in time.
- Argument 'c' is the velocity of the wave
- Optional argument 'sav': if a string is given, then saves data and figures.

"""
function figure8(c::Real;sav=[])
	param = ( μ  = 1, ϵ  = 1, c=c ,
	    L  = 10*π, N  = 2^10,
		T  = 1, dt = 1/2000, ns=1
		)
	mesh = Mesh(param)
	function solη(x,c,ϵ,μ)
		(c^2-1)/ϵ*sech(sqrt(3*(c^2-1)/(c^2)/μ)/2*x)^2
	end
	η(c,x0)= solη.(mesh.x .-x0,c,param.ϵ,param.μ)
	u(c,x0)= c*η(c,x0)./(1 .+ param.ϵ*η(c,x0))

	(ηinit,uinit) = SolitaryWaveWhithamGreenNaghdi(
			mesh, param, η(c,0);
			method=2, α = 1,
			tol =  1e-13, max_iter=10,
			iterative = false, q=1,
			verbose = true, SGN = false)
	(ηfinexact,ufinexact) = SolitaryWaveWhithamGreenNaghdi(
			mesh, param, η(c,c*param.T);
			method=2, α = 1,
			tol =  1e-13, max_iter=10,
			iterative = false, q=1,
			verbose = true, SGN = false)

	hinit = 1 .+ param.ϵ*ηinit
	F₁ 	= tanh.(sqrt(param.μ)*abs.(mesh.k))./(sqrt(param.μ)*abs.(mesh.k))
	F₁[1]= 1                 # Differentiation
	F₀   = 1im * sqrt.(3*(1 ./F₁ .- 1)).*sign.(mesh.k)
	DxF(v) = real.(ifft(F₀ .* fft(v)))
	vinit= uinit - 1/3 ./hinit .* (DxF(hinit.^3 .*DxF(uinit)))

	init     = Init(mesh,ηinit,vinit)
	model = WhithamGreenNaghdiKlein(param;SGN=false, ktol=0*1e-16, iterate=true, precond = false)
	problem = Problem(model, init, param)

	@time solve!( problem )

	E(η,u,v) = sum(η.^2 .+ (1 .+ param.ϵ*η).*u.*v)
	E1(η,u) = sum(η.^2 .+ (1 .+ param.ϵ*η).*u.^2 .+1/3*(1 .+ param.ϵ*η).^3 .*(DxF(u).^2))
	E2(η,u) = sum(η.^2 .+ (1 .+ param.ϵ*η).*u.^2 .-1/3*u.*DxF((1 .+ param.ϵ*η).^3 .*(DxF(u))))
	dE(η1,u1,v1,η2,u2,v2) = sum(η1.^2-η2.^2) + sum((1 .+ param.ϵ*η1).*u1.*v1 - (1 .+ param.ϵ*η2).*u2.*v2)
	dE1(η1,u1,η2,u2) = sum(η1.^2-η2.^2) + sum((1 .+ param.ϵ*η1).*u1.^2 - (1 .+ param.ϵ*η2).*u2.^2) -
			1/3*sum(u1.*DxF((1 .+ param.ϵ*η1).^3 .*(DxF(u1)))-u2.*DxF((1 .+ param.ϵ*η2).^3 .*(DxF(u2))))
	dE2(η1,u1,η2,u2) = sum(η1.^2-η2.^2) + sum((1 .+ param.ϵ*η1).*u1.^2 - (1 .+ param.ϵ*η2).*u2.^2) +
			1/3*sum((1 .+ param.ϵ*η1).^3 .*(DxF(u1).^2) - (1 .+ param.ϵ*η2).^3 .*(DxF(u2).^2))

	(ηfin,vfin,ufin)   =  mapfro(model,last(problem.data.U))
	print(string("normalized error: ",dE(ηinit,uinit,vinit,ηfin,ufin,vfin)/E(ηinit,uinit,vinit),"\n"))
	print(string("normalized error: ",dE1(ηinit,uinit,ηfin,ufin)/E1(ηinit,uinit),"\n"))
	print(string("normalized error: ",dE2(ηinit,uinit,ηfin,ufin)/E2(ηinit,uinit),"\n"))

	(ηinit,vinit,uinit) = mapfro(model,first(problem.data.U))
	print(string("normalized error: ",dE(ηinit,uinit,vinit,ηfin,ufin,vfin)/E(ηinit,uinit,vinit),"\n"))
	print(string("normalized error: ",dE1(ηinit,uinit,ηfin,ufin)/E1(ηinit,uinit),"\n"))
	print(string("normalized error: ",dE2(ηinit,uinit,ηfin,ufin)/E2(ηinit,uinit),"\n"))

	plt = plot(layout=(1,2))
	plot!(plt[1,1],fftshift(mesh.k),fftshift(log10.(abs.(fft(vinit))));
			title = "frequency",
			label = "initial")
	plot!(plt[1,1],fftshift(mesh.k),fftshift(log10.(abs.(fft(vfin))));
			label = "final")
	plot!(plt[1,2],mesh.x,ufin-ufinexact;
			title = string("error at time t=",problem.times.tfin),
			label="difference")
	display(plt)

	if sav != []
		savefig(string("fig8-",sav,".pdf"));
		#p = plot(layout=(2,1))
		#fig_problem!(p,problem)
		#savefig(string("fig8a-",sav,".pdf"));
		save(problem,string("fig8-",sav,".pdf"));
	end
end
#figure8(c) # c is the velocity


#------ Figures 9 to 11

"""
	figure10(p,c;sav)

perturbed SGN solitary wave integrated in time.
- Argument 'p' (1,2,3 or 4) is the type of perturbation
- Argument 'c' is the velocity of the wave
- Optional argument 'sav': if a string is given, then saves data and figures.

"""
p=4
c=4
sav="c4tp01"
function figure10(p::Int,c::Real;sav=[])
	if p == 1
		λ = 0.99
	elseif p == 2
		λ = 1.01
	elseif p == 3
		λ = -0.01
	elseif p == 4
		λ = +0.01
	end
	param = ( μ = 1, ϵ = 1, c = c, λ = λ,
				N  = 2^12,
	            L  = 30*π,
	            T  = 20,
	            dt = 20/10^4, ns=100,
				dealias = 0,
				SGN=true,
				ktol=0*1e-6,
				gtol=1e-12,
				iterate=true,
				precond = false)
	function solη(x,c,ϵ,μ)
		(c^2-1)/ϵ*sech(sqrt(3*(c^2-1)/(c^2)/μ)/2*x)^2
	end
	function solu(x,c,ϵ,μ)
		c*solη(x,c,ϵ,μ) / (1 + ϵ * solη(x,c,ϵ,μ))
	end
	mesh=Mesh(param)
	ηGN= solη.(mesh.x,c,param.ϵ,param.μ)
	hGN=(1 .+ param.ϵ*ηGN)
	uGN=solu.(mesh.x,c,param.ϵ,param.μ)
	if p == 1 || p == 2
		uGN= λ*uGN
	elseif p == 3 || p == 4
		uGN .+= λ*exp.(-mesh.x.^2)
	end
	Dx(v) = real.(ifft(1im*mesh.k .* fft(v)))
	vGN= uGN - param.μ/3 ./hGN .* (Dx(hGN.^3 .*Dx(uGN)))

	init = Init(mesh,ηGN,vGN)
	model = WhithamGreenNaghdi(param;SGN=param.SGN, dealias = param.dealias, ktol=param.ktol, gtol=param.gtol, iterate=param.iterate, precond = param.precond)
	problem = Problem(model, init, param)

	@time solve!( problem )

	if sav != [] save(problem,string("fig10-",sav)); end

	(ηfin,vfin,ufin)   =  mapfro(model,last(problem.data.U))
	(ηinit,vinit,uinit) = mapfro(model,first(problem.data.U))


	E(η) = sum(η.^2 ./ (1 .+ η) .+ param.μ/3*(1 .+ param.ϵ*η).^3 .*(Dx(η ./ (1 .+ η)).^2))
	E(η,u,v) = sum(η.^2 .+ (1 .+ param.ϵ*η).*u.*v)
	E(η,u) = sum(η.^2 .+ (1 .+ param.ϵ*η).*u.^2 .+param.μ/3*(1 .+ param.ϵ*η).^3 .*(Dx(u).^2))
	dE(η1,u1,v1,η2,u2,v2) = sum(η1.^2-η2.^2) + sum((1 .+ param.ϵ*η1).*u1.*v1 - (1 .+ param.ϵ*η2).*u2.*v2)
	dE(η1,u1,η2,u2) = sum(η1.^2-η2.^2) + sum((1 .+ param.ϵ*η1).*u1.^2 - (1 .+ param.ϵ*η2).*u2.^2) -
			param.μ/3*sum(u1.*Dx((1 .+ param.ϵ*η1).^3 .*(Dx(u1)))-u2.*Dx((1 .+ param.ϵ*η2).^3 .*(Dx(u2))))

	print(string("normalized error: ",dE(ηGN,uGN,vGN,ηfin,ufin,vfin)/E(ηGN,uGN,vGN),"\n"))
	print(string("normalized error: ",dE(ηGN,uGN,ηfin,ufin)/E(ηGN,uGN),"\n"))
	print(string("normalized error: ",dE(ηinit,uinit,vinit,ηfin,ufin,vfin)/E(ηinit,uinit,vinit),"\n"))
	print(string("normalized error: ",dE(ηinit,uinit,ηfin,ufin)/E(ηinit,uinit),"\n"))

	plt = plot(layout=(1,2))
	plot!(plt[1,1],fftshift(mesh.k),fftshift(log10.(abs.(fft(ηinit))));
			title = "frequency",
			label = "initial")
	plot!(plt[1,1],fftshift(mesh.k),fftshift(log10.(abs.(fft(ηfin))));
			label = "final")
	plot!(plt[1,2],mesh.x,[ηfin solη.(mesh.x.-c*param.T,c,param.ϵ,param.μ)];
			title = string("at time t=",problem.times.tfin),
			label=["zeta" "unperturbed soliton"])
	plot!(plt[1,2],mesh.x,[ufin solu.(mesh.x.-c*param.T,c,param.ϵ,param.μ)];
			title = string("at time t=",problem.times.tfin),
			label=["u" "unperturbed soliton"])
	display(plt)
	if sav != [] savefig(string("fig9a-",sav,".pdf")); end

	ts = problem.times.ts
	x = mesh.x
	k = mesh.k
	us=zeros(param.N,length(ts));
	h = similar(Complex.(x))
	fftu= Complex.(h)
	Precond = Diagonal( 1 .+ param.μ*k.^2 )

	function compute(fftη,fftv)
		h .= 1 .+ ifft(fftη)
		function LL(hatu)
			hatu- param.μ/3 *fft( 1 ./h .* ifft( 1im*k .* fft( h.^3 .* ifft( 1im*k .* hatu ) ) ) )
		end
		fftu.= gmres( LinearMap(LL, length(h); issymmetric=false, ismutating=false) , fftv ;
					Pl = Precond )
		return real.(ifft(fftu))
	end
	@showprogress 1 for i in 1:length(ts)
		#us[:,i].=compute(problem.data.U[i][:,1],problem.data.U[i][:,2])
		us[:,i].=real.(ifft(problem.data.U[i][:,1]))
	end

	# plt=plot()
	# my_cg = cgrad([:blue,:green])
	# surface!(plt,x,ts,us',view_angle=(20,30), color = my_cg)
	# display(plt)
	# if sav != [] savefig(string("fig9b-",sav,".pdf")); end

	plt = plot()
	plot!(plt,ts,maximum(abs.(us),dims=1)',
			title="L infty norm",
			label="")
	display(plt)

	if sav != [] savefig(string("fig10-",sav,".pdf")); end
end
#figure10(4,10,sav="test1")  # use with p = 1,2,3 or 4 (type of perturbation); and c=2 (velocity)
#figure12(4,10)  # use with p = 1,2,3 or 4 (type of perturbation); and c=2 (velocity)

#------ Figure 12

"""
	figure12(p,c;sav)

perturbed WGN solitary wave integrated in time.
- Argument 'p' (1,2,3 or 4) is the type of perturbation
- Argument 'c' is the velocity of the wave
- Optional argument 'sav': if a string is given, then saves data and figures.

"""
function figure12(p::Int,c::Real;sav=[])
	if p == 1
		λ = 0.99
	elseif p == 2
		λ = 1.01
	end
	param = ( μ  = 1, ϵ  = 1, c=c,
				N  = 2^10,
	            L  = 10*π,
	            T  = 1.5,
	            dt = 1.5/10^4, ns=100)
	mesh=Mesh(param)

	function solη(x,c,ϵ,μ)
		(c^2-1)/ϵ*sech(sqrt(3*(c^2-1)/(c^2)/μ)/2*x)^2
	end
	η(c,x0)= solη.(mesh.x .-x0,c,param.ϵ,param.μ)
	u(c,x0)= c*η(c,x0)./(1 .+ param.ϵ*η(c,x0))

	(ηinit,uinit) = SolitaryWaveWhithamGreenNaghdi(
			mesh, param, η(c,0);
			method=3, α = 1,
			tol =  1e-13, max_iter=10,
			iterative = false, q=1,
			verbose = true, SGN = false)

	hinit = 1 .+ param.ϵ*ηinit

	if p == 1 || p == 2
		uinit= λ*uinit
	elseif p == 3
		uinit.-= 0.01*exp.(-mesh.x.^2)
	elseif p == 4
		uinit.+= exp.(-mesh.x.^2)
	end

	F₁ 	= tanh.(sqrt(param.μ)*abs.(mesh.k))./(sqrt(param.μ)*abs.(mesh.k))
	F₁[1]= 1                 # Differentiation
	F₀   = 1im * sqrt.(3*(1 ./F₁ .- 1)).*sign.(mesh.k)
	DxF(v) = real.(ifft(F₀ .* fft(v)))
	vinit= uinit - 1/3 ./hinit .* (DxF(hinit.^3 .*DxF(uinit)))


	init     = Init(mesh,ηinit,vinit)
	model = WhithamGreenNaghdi(param;SGN=false, ktol=0, iterate=false, precond = false)
	problem = Problem(model, init, param)

	@time solve!( problem )
	if sav != []
		save(problem,string("fig12-",sav))
	end

	(ηfin,vfin,ufin)   =  mapfro(model,last(problem.data.U))
	(ηini,vini,uini) = mapfro(model,first(problem.data.U))


	E(η,u,v) = sum(η.^2 .+ (1 .+ param.ϵ*η).*u.*v)
	E(η,u) = sum(η.^2 .+ (1 .+ param.ϵ*η).*u.^2 .+param.μ/3*(1 .+ param.ϵ*η).^3 .*(DxF(u).^2))
	dE(η1,u1,v1,η2,u2,v2) = sum(η1.^2-η2.^2) + sum((1 .+ param.ϵ*η1).*u1.*v1 - (1 .+ param.ϵ*η2).*u2.*v2)
	dE(η1,u1,η2,u2) = sum(η1.^2-η2.^2) + sum((1 .+ param.ϵ*η1).*u1.^2 - (1 .+ param.ϵ*η2).*u2.^2) -
			param.μ/3*sum(u1.*DxF((1 .+ param.ϵ*η1).^3 .*(DxF(u1)))-u2.*DxF((1 .+ param.ϵ*η2).^3 .*(DxF(u2))))


	print(string("normalized error: ",dE(ηinit,uinit,vinit,ηfin,ufin,vfin)/E(ηinit,uinit,vinit),"\n"))
	print(string("normalized error: ",dE(ηinit,uinit,ηfin,ufin)/E(ηinit,uinit),"\n"))
	print(string("normalized error: ",dE(ηini,uini,vini,ηfin,ufin,vfin)/E(ηini,uini,vini),"\n"))
	print(string("normalized error: ",dE(ηini,uini,ηfin,ufin)/E(ηini,uini),"\n"))


	plt = plot(layout=(1,2))
	plot!(plt[1,1],fftshift(mesh.k),fftshift(log10.(abs.(fft(vini))));
			title = "frequency",
			label = "initial")
	plot!(plt[1,1],fftshift(mesh.k),fftshift(log10.(abs.(fft(vfin))));
			label = "final")
	plot!(plt[1,2],mesh.x,ufin;
			title = string("at time t=",problem.times.tfin),
			label=["u"])
	display(plt)

	if sav != [] savefig(string("fig12a-",sav,".pdf")); 	end


	ts = problem.times.ts
	x = mesh.x
	k = mesh.k
	us=zeros(param.N,length(ts));
	h = similar(Complex.(x))
	fftu= Complex.(h)
	Precond = Diagonal( 1 .+ param.μ*k.^2 )

	function compute(fftη,fftv)
		h .= 1 .+ ifft(fftη)
		function LL(hatu)
			hatu- 1/3 *fft( 1 ./h .* ifft( F₀ .* fft( h.^3 .* ifft( F₀ .* hatu ) ) ) )
		end
		fftu.= gmres( LinearMap(LL, length(h); issymmetric=false, ismutating=false) , fftv ;
					Pl = Precond )
		return real.(ifft(fftu))
	end
	@showprogress 1 for i in 1:length(ts)
		us[:,i].=compute(problem.data.U[i][:,1],problem.data.U[i][:,2])
	end

	if sav != []
		plt=plot()
		my_cg = cgrad([:blue,:green])
		surface!(plt,x,ts,us',view_angle=(20,30), color = my_cg)
		display(plt)
		savefig(string("fig12b-",sav,".pdf"));
	end
	plt = plot()
	plot!(plt,ts,maximum(abs.(us),dims=1)',
			title="L infty norm",
			label="u")
	#plot!(plt[1,2],ts,[normalize(mini,Inf) normalize(L2,Inf)],
	#		title=["minimizer functional" "L2 norm"],
	#		label="")
	display(plt)

	if sav != [] savefig(string("fig12c-",sav,".pdf")); end
end
#figure12(p,c)  # use with p = 1,2,3 or 4 (type of perturbation); and c=2 (velocity)
