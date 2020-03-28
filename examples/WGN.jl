"""
Reproduces some figures in the work of C. Klein and V. Duchêne
"""
#using ShallowWaterModels
include("../src/dependencies.jl")


#---- Figure 7
function figure7()
	param = ( μ  = 1, ϵ  = 1, c=2,
				N  = 2^9,
	            L  = 10*π,
	            T  = 1,
	            dt = 1/2000, ns=500)
	function solη(x,c,ϵ,μ)
		(c^2-1)/ϵ*sech(sqrt(3*(c^2-1)/(c^2)/μ)/2*x)^2
	end
	function solu(x,c,ϵ,μ)
		c*solη(x,c,ϵ,μ) / (1 + ϵ * solη(x,c,ϵ,μ))
	end
	mesh=Mesh(param)
	ηGN= solη.(mesh.x,param.c,param.ϵ,param.μ)
	hGN=(1 .+ param.ϵ*ηGN)
	uGN= solu.(mesh.x,param.c,param.ϵ,param.μ)
	Dx(v) = real.(ifft(1im*mesh.k .* fft(v)))
	vGN= uGN - param.μ/3 ./hGN .* (Dx(hGN.^3 .*Dx(uGN)))

	init     = Init(mesh,ηGN,vGN)
	model = WhithamGreenNaghdi(param;SGN=true, ktol=0, iterate=true, precond = false)
	problem = Problem(model, init, param)

	@time solve!( problem )

	(ηfin,vfin,ufin)   =  mapfro(model,last(problem.data.U))
	(ηinit,vinit,uinit) = mapfro(model,first(problem.data.U))


	E(η,v,u) = sum(η.^2 .+ (1 .+ param.ϵ*η).*u.*v)
	E1(η,u) = sum(η.^2 .+ (1 .+ param.ϵ*η).*u.^2 .+param.μ/3*(1 .+ param.ϵ*η).^3 .*(Dx(u).^2))
	E2(η,u) = sum(η.^2 .+ (1 .+ param.ϵ*η).*u.^2 .-param.μ/3*u.*DxF((1 .+ param.ϵ*η).^3 .*(Dx(u))))
	dE(η1,u1,v1,η2,u2,v2) = sum(η1.^2-η2.^2) + sum((1 .+ param.ϵ*η1).*u1.*v1 - (1 .+ param.ϵ*η2).*u2.*v2)
	dE1(η1,u1,η2,u2) = sum(η1.^2-η2.^2) + sum((1 .+ param.ϵ*η1).*u1.^2 - (1 .+ param.ϵ*η2).*u2.^2) +
			param.μ/3*sum(u1.*Dx((1 .+ param.ϵ*η1).^3 .*(Dx(u1)))-u2.*Dx((1 .+ param.ϵ*η2).^3 .*(Dx(u2))))
	dE2(η1,u1,η2,u2) = sum(η1.^2-η2.^2) + sum((1 .+ param.ϵ*η1).*u1.^2 - (1 .+ param.ϵ*η2).*u2.^2) +
			param.μ/3*sum((1 .+ param.ϵ*η1).^3 .*(Dx(u1).^2) - (1 .+ param.ϵ*η2).^3 .*(Dx(u2).^2))

	print(string("final energy minus initial energy: ",dE(ηGN,uGN,vGN,ηfin,ufin,vfin),"\n"))
	print(string("final energy minus initial energy: ",dE1(ηGN,uGN,ηfin,ufin),"\n"))
	print(string("final energy minus initial energy: ",dE2(ηGN,uGN,ηfin,ufin),"\n"))

	print(string("final energy minus initial energy: ",dE(ηinit,uinit,vinit,ηfin,ufin,vfin),"\n"))
	print(string("final energy minus initial energy: ",dE1(ηinit,uinit,ηfin,ufin),"\n"))
	print(string("final energy minus initial energy: ",dE2(ηinit,uinit,ηfin,ufin),"\n"))


	#p = plot(layout=(2,1))
	#fig_problem!( p,problem)

	plt = plot(layout=(1,2))
	plot!(plt[1,1],fftshift(mesh.k),fftshift(log10.(abs.(fft(vinit))));
			title = "frequency",
			label = "initial")
	plot!(plt[1,1],fftshift(mesh.k),fftshift(log10.(abs.(fft(vfin))));
			label = "final")
	plot!(plt[1,2],mesh.x,ufin-solu.(mesh.x.-param.c*param.T,param.c,param.ϵ,param.μ);
			title = string("error at time t=",problem.times.tfin),
			label="difference")
	display(plt)

	savefig("fig7iv.pdf");
	#save(problem,"fig7")
end
figure7()
#------ Figure 8
function figure8()
	param = ( μ  = 1, ϵ  = 1, c = 2,
	    L  = 10*π, N  = 2^10,
		T  = 0.1, dt = 1/2000, ns=50
		)
	mesh = Mesh(param)
	function solη(x,c,ϵ,μ)
		(c^2-1)/ϵ*sech(sqrt(3*(c^2-1)/(c^2)/μ)/2*x)^2
	end
	η(c,x0)= solη.(mesh.x .-x0,c,param.ϵ,param.μ)
	u(c,x0)= c*η(c,x0)./(1 .+ param.ϵ*η(c,x0))

	(ηinit,uinit) = SolitaryWaveWhithamGreenNaghdi(
			mesh, param, η(param.c,0);
			method=2, α = 1,
			tol =  1e-13, max_iter=10,
			iterative = true, q=1,
			verbose = false, SGN = false)
	(ηfinexact,ufinexact) = SolitaryWaveWhithamGreenNaghdi(
			mesh, param, η(param.c,param.c*param.T);
			method=2, α = 1,
			tol =  1e-13, max_iter=10,
			iterative = true, q=1,
			verbose = false, SGN = false)

	hinit = 1 .+ param.ϵ*ηinit
	F₁ 	= tanh.(sqrt(param.μ)*abs.(mesh.k))./(sqrt(param.μ)*abs.(mesh.k))
	F₁[1]= 1                 # Differentiation
	F₀   = 1im * sqrt.(3*(1 ./F₁ .- 1)).*sign.(mesh.k)
	DxF(v) = real.(ifft(F₀ .* fft(v)))
	vinit= uinit - 1/3 ./hinit .* (DxF(hinit.^3 .*DxF(uinit)))

	init     = Init(mesh,ηinit,vinit)
	model = WhithamGreenNaghdiKlein(param;SGN=false, ktol=0*1e-16, iterate=true, precond = false)
	problem = Problem(model2, init, param)

	@time solve!( problem )

	E(η,v,u) = sum(η.^2 .+ (1 .+ param.ϵ*η).*u.*v)
	E1(η,u) = sum(η.^2 .+ (1 .+ param.ϵ*η).*u.^2 .+1/3*(1 .+ param.ϵ*η).^3 .*(DxF(u).^2))
	E2(η,u) = sum(η.^2 .+ (1 .+ param.ϵ*η).*u.^2 .-1/3*u.*DxF((1 .+ param.ϵ*η).^3 .*(DxF(u))))
	dE(η1,u1,v1,η2,u2,v2) = sum(η1.^2-η2.^2) + sum((1 .+ param.ϵ*η1).*u1.*v1 - (1 .+ param.ϵ*η2).*u2.*v2)
	dE1(η1,u1,η2,u2) = sum(η1.^2-η2.^2) + sum((1 .+ param.ϵ*η1).*u1.^2 - (1 .+ param.ϵ*η2).*u2.^2) +
			1/3*sum(u1.*DxF((1 .+ param.ϵ*η1).^3 .*(DxF(u1)))-u2.*DxF((1 .+ param.ϵ*η2).^3 .*(DxF(u2))))
	dE2(η1,u1,η2,u2) = sum(η1.^2-η2.^2) + sum((1 .+ param.ϵ*η1).*u1.^2 - (1 .+ param.ϵ*η2).*u2.^2) +
			1/3*sum((1 .+ param.ϵ*η1).^3 .*(DxF(u1).^2) - (1 .+ param.ϵ*η2).^3 .*(DxF(u2).^2))

	(ηfin,vfin,ufin)   =  mapfro(model,last(problem.data.U))
	print(string("final energy minus initial energy: ",dE(ηinit,uinit,vinit,ηfin,ufin,vfin),"\n"))
	print(string("final energy minus initial energy: ",dE1(ηinit,uinit,ηfin,ufin),"\n"))
	print(string("final energy minus initial energy: ",dE2(ηinit,uinit,ηfin,ufin),"\n"))

	(ηinit,vinit,uinit) = mapfro(model,first(problem.data.U))
	print(string("final energy minus initial energy: ",dE(ηinit,uinit,vinit,ηfin,ufin,vfin),"\n"))
	print(string("final energy minus initial energy: ",dE1(ηinit,uinit,ηfin,ufin),"\n"))
	print(string("final energy minus initial energy: ",dE2(ηinit,uinit,ηfin,ufin),"\n"))

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

	#savefig("fig8s.pdf");
	#save(problem,"fig8s")
end

figure8()
