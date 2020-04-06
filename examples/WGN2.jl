export figure7,figure8,figure10,figure12

"""
Reproduces some figures in the work of C. Klein and V. Duchêne
"""
#using ShallowWaterModels
include("../src/dependencies.jl")


#---- Figure 14
"""
	figure14(ε,SGN;sav)

Dispersive shock wave for SGN or WGN.
- Argument `ε` is the square root of the shallowness parameter
- Solves SGN if `SGN` is `true`, WGN otherwise.
- Optional argument 'sav': if a string is given, then saves data and figures.

"""
function figure14(ε::Real,SGN::Bool;sav=[])
	param = ( μ  = ε^2, ϵ  = 1,
				N  = 2^11,
	            L  = 3*π,
	            T  = 5,
	            dt = 5/10000, ns=1000)

	mesh=Mesh(param)
	η= exp.(-(mesh.x .+3).^2)
	h=1 .+ param.ϵ*η
	u= 2*sqrt.(h) .-2
	if SGN == true
		F₀ = sqrt(param.μ)*1im*mesh.k
	else
		F₁ 	= tanh.(sqrt(param.μ)*abs.(mesh.k))./(sqrt(param.μ)*abs.(mesh.k))
		F₁[1]= 1                 # Differentiation
		F₀   = 1im * sqrt.(3*(1 ./F₁ .- 1)).*sign.(mesh.k)
	end
	DxF(v) = real.(ifft(F₀ .* fft(v)))
	v= u - 1/3 ./h .* (DxF(h.^3 .*DxF(u)))
	v= 2*sqrt.(h) .-2

	init     = Init(mesh,η,v)
	model = WhithamGreenNaghdi(param;SGN=SGN, ktol=0, gtol=1e-12, iterate=true, precond = false)
	problem = Problem(model, init, param)

	@time solve!( problem )
	#
	# (ηfin,vfin,ufin)   =  mapfro(model,last(problem.data.U))
	# (ηinit,vinit,uinit) = mapfro(model,first(problem.data.U))
	#
	#
	# E(η,u,v) = sum(η.^2 .+ (1 .+ param.ϵ*η).*u.*v)
	# E1(η,u) = sum(η.^2 .+ (1 .+ param.ϵ*η).*u.^2 .+1/3*(1 .+ param.ϵ*η).^3 .*(DxF(u).^2))
	# E2(η,u) = sum(η.^2 .+ (1 .+ param.ϵ*η).*u.^2 .-1/3*u.*DxF((1 .+ param.ϵ*η).^3 .*(DxF(u))))
	# dE(η1,u1,v1,η2,u2,v2) = sum(η1.^2-η2.^2) + sum((1 .+ param.ϵ*η1).*u1.*v1 - (1 .+ param.ϵ*η2).*u2.*v2)
	# dE1(η1,u1,η2,u2) = sum(η1.^2-η2.^2) + sum((1 .+ param.ϵ*η1).*u1.^2 - (1 .+ param.ϵ*η2).*u2.^2) -
	# 		1/3*sum(u1.*DxF((1 .+ param.ϵ*η1).^3 .*(DxF(u1)))-u2.*DxF((1 .+ param.ϵ*η2).^3 .*(DxF(u2))))
	# dE2(η1,u1,η2,u2) = sum(η1.^2-η2.^2) + sum((1 .+ param.ϵ*η1).*u1.^2 - (1 .+ param.ϵ*η2).*u2.^2) +
	# 		1/3*sum((1 .+ param.ϵ*η1).^3 .*(DxF(u1).^2) - (1 .+ param.ϵ*η2).^3 .*(DxF(u2).^2))
	#
	# print(string("normalized error: ",dE(η,u,v,ηfin,ufin,vfin)/E(η,u,v),"\n"))
	# print(string("normalized error: ",dE1(η,u,ηfin,ufin)/E1(η,u),"\n"))
	# # print(string("normalized error: ",dE2(η,u,ηfin,ufin)/E2(η,u),"\n"))
	#
	# print(string("normalized error: ",dE(ηinit,uinit,vinit,ηfin,ufin,vfin)/E(ηinit,uinit,vinit),"\n"))
	# print(string("normalized error: ",dE1(ηinit,uinit,ηfin,ufin)/E1(ηinit,uinit),"\n"))
	# print(string("normalized error: ",dE2(ηinit,uinit,ηfin,ufin)/E2(ηinit,uinit),"\n"))

	fftηfin=last(problem.data.U)[:,1]

	plt = plot(layout=(1,2))
	#plot!(plt[1,1],fftshift(mesh.k),fftshift(log10.(abs.(fft(ηinit))));
	#		title = "frequency",
	#		label = "initial")
	plot!(plt[1,1],fftshift(mesh.k),fftshift(log10.(abs.(fftηfin)));
			title = "frequency",
			label = "final")
	#plot!(plt[1,1],fftshift(mesh.k),fftshift(log10.(abs.(fftηfin)));
	#		label = "final 2")
	plot!(plt[1,2],mesh.x,real.(ifft(fftηfin));
			title = string("surface time t=",problem.times.tfin),
			label="eta")
	xlims!(plt[1,2],(0,8))
	display(plt)

	if sav != []
		savefig(string("fig14-",sav,".pdf"));
		#p = plot(layout=(2,1))
		#fig_problem!(p,problem)
		#savefig(string("fig7a-",sav,".pdf"));
		save(problem,string("fig14-",sav));
	end
end

#figure14(0.01,true;sav="SGN")

#---- Figure 14
"""
	figure14(ε,SGN;sav)

Dispersive shock wave for SGN or WGN.
- Argument `ε` is the square root of the shallowness parameter
- Solves SGN if `SGN` is `true`, WGN otherwise.
- Optional argument 'sav': if a string is given, then saves data and figures.

"""
function figure17(ε::Real,SGN::Bool;sav=[])
	param = ( μ  = ε^2, ϵ  = 1,
				N  = 2^11,
	            L  = 5*π,
	            T  = 10,
	            dt = 10/10^3, ns=10,
				dealias = 1,
				SGN=SGN,
				ktol=0*1e-6,
				gtol=1e-12,
				iterate=true, precond = false)

	mesh=Mesh(param)
	function krasny!(v)
		v[abs.(v).<1e-14].=0
		return v
	end
	Dx(v) = real.(ifft( 1im*mesh.k.* krasny!(fft(v))))
	Dx2(v) = real.(ifft( -(mesh.k.^2) .* fft(v)))
	Dx2(v) = Dx(Dx(v))
	ϵ = param.ϵ;μ=param.μ;
	w = - (mesh.x).* exp.(-(mesh.x).^2)
	u = w .+ μ/12 *Dx2(w) .+ μ*ϵ/6* w.*Dx2(w)
	η = u .+ ϵ/4*u.^2 .- μ/6* Dx2( u+3*ϵ/4* u.^2) .- μ*ϵ/6 * u .* Dx2(u) .- 5*μ*ϵ/48 * Dx(u).^2
	h=1 .+ param.ϵ*η
	if SGN == true
		F₀ = sqrt(param.μ)*1im*mesh.k
	else
		F₁ 	= tanh.(sqrt(param.μ)*abs.(mesh.k))./(sqrt(param.μ)*abs.(mesh.k))
		F₁[1]= 1                 # Differentiation
		F₀   = 1im * sqrt.(3*(1 ./F₁ .- 1)).*sign.(mesh.k)
	end
	DxF(v) = real.(ifft(F₀ .* fft(v)))
	v= u - 1/3 ./h .* (DxF(h.^3 .*DxF(u)))

	init     = Init(mesh,η,v)
	model = WhithamGreenNaghdi(param;SGN=param.SGN, dealias = param.dealias, ktol=param.ktol, gtol=param.gtol, iterate=param.iterate, precond = param.precond)
	problem = Problem(model, init, param)

	@time solve!( problem )

	if sav != [] save(problem,string("fig17-",sav)); end

	fftηfin=last(problem.data.U)[:,1]

	plt = plot(layout=(1,2))
	#plot!(plt[1,1],fftshift(mesh.k),fftshift(log10.(abs.(fft(ηinit))));
	#		title = "frequency",
	#		label = "initial")
	plot!(plt[1,1],fftshift(mesh.k),fftshift(log10.(abs.(fftηfin)));
			title = "frequency",
			label = "final")
	#plot!(plt[1,1],fftshift(mesh.k),fftshift(log10.(abs.(fftηfin)));
	#		label = "final 2")
	plot!(plt[1,2],mesh.x,real.(ifft(fftηfin));
			title = string("surface time t=",problem.times.tfin),
			label="eta")
	display(plt)


	ts = problem.times.ts
	us=zeros(param.N,length(ts));
	@showprogress 1 for i in 1:length(ts)
		#us[:,i].=compute(problem.data.U[i][:,1],problem.data.U[i][:,2])
		us[:,i].=real.(ifft(problem.data.U[i][:,1]))
	end

	plt=plot()
	my_cg = cgrad([:blue,:green])
	surface!(plt,mesh.x,ts,us',view_angle=(20,30), color = my_cg)
	display(plt)
	if sav != [] savefig(string("fig17-",sav,".pdf")); end

	# plt = plot()
	# plot!(plt,ts,maximum(abs.(us),dims=1)',
	# 		title="L infty norm",
	# 		label="")
	# display(plt)


	anim = @animate for i in range(1,length(ts))
		plt = plot(layout=(1,2))
		plot!(plt[1,1],fftshift(mesh.k),fftshift(log10.(abs.(problem.data.U[i][:,1])));
				title = "frequency",
				label = "")
		plot!(plt[1,2],mesh.x,us[:,i];
				title = string("surface at time t=",problem.times.ts[i]),
				label="")
		#display(plt)
		#ylims!(plt[1,2],-1,18)
		#xlims!(plt[1,2],-10,75)
		#next!(prog)
	end

	gif(anim, string("anim-",sav,".gif"), fps=15); nothing

	# if sav != []
	#
	# 	savefig(string("fig17-",sav,".pdf"));
	# 	p = plot(layout=(2,1))
	# 	fig_problem!(p,problem)
	# 	savefig(string("fig18-",sav,".pdf"));
	# 	save(problem,string("fig17-",sav));
	# end
end

function figure19(λ::Real,SGN::Bool;sav=[])
	param = ( μ  = 1, ϵ  = 1,
				N  = 2^12,
	            L  = 3*π,
	            T  = 0.6,
	            dt = 0.6/10^4, ns=100,
				dealias = 1,
				SGN=SGN,
				ktol=0*1e-6,
				gtol=1e-12,
				iterate=true, precond = true)

	mesh=Mesh(param)
	η = -λ*exp.(-mesh.x.^2)
	u = -mesh.x.*exp.(-mesh.x.^2)
	h=1 .+ param.ϵ*η
	if SGN == true
		F₀ = sqrt(param.μ)*1im*mesh.k
	else
		F₁ 	= tanh.(sqrt(param.μ)*abs.(mesh.k))./(sqrt(param.μ)*abs.(mesh.k))
		F₁[1]= 1                 # Differentiation
		F₀   = 1im * sqrt.(3*(1 ./F₁ .- 1)).*sign.(mesh.k)
	end
	DxF(v) = real.(ifft(F₀ .* fft(v)))
	v= u - 1/3 ./h .* (DxF(h.^3 .*DxF(u)))

	init     = Init(mesh,η,v)
	model = WhithamGreenNaghdi(param;SGN=param.SGN, dealias = param.dealias, ktol=param.ktol, gtol=param.gtol, iterate=param.iterate, precond = param.precond)
	problem = Problem(model, init, param)

	@time solve!( problem )

	if sav != [] save(problem,string("fig19-",sav)); end

	fftηfin=last(problem.data.U)[:,1]

	plt = plot(layout=(1,2))
	#plot!(plt[1,1],fftshift(mesh.k),fftshift(log10.(abs.(fft(ηinit))));
	#		title = "frequency",
	#		label = "initial")
	plot!(plt[1,1],fftshift(mesh.k),fftshift(log10.(abs.(fftηfin)));
			title = "frequency",
			label = "final")
	#plot!(plt[1,1],fftshift(mesh.k),fftshift(log10.(abs.(fftηfin)));
	#		label = "final 2")
	plot!(plt[1,2],mesh.x,real.(ifft(fftηfin));
			title = string("surface time t=",problem.times.tfin),
			label="eta")
	display(plt)


	ts = problem.times.ts
	us=zeros(param.N,length(ts));
	@showprogress 1 for i in 1:length(ts)
		#us[:,i].=compute(problem.data.U[i][:,1],problem.data.U[i][:,2])
		us[:,i].=real.(ifft(problem.data.U[i][:,1]))
	end

	plt=plot()
	my_cg = cgrad([:blue,:green])
	surface!(plt,mesh.x,ts,us',view_angle=(20,30), color = my_cg)
	display(plt)
	if sav != [] savefig(string("fig19t-",sav,".pdf")); end

	# plt = plot()
	# plot!(plt,ts,maximum(abs.(us),dims=1)',
	# 		title="L infty norm",
	# 		label="")
	# display(plt)
	function plotfig(i)
		plt = plot(layout=(1,2))
		plot!(plt[1,1],fftshift(mesh.k),fftshift(log10.(abs.(problem.data.U[i][:,1])));
				title = "frequency",
				label = "")
		plot!(plt[1,2],mesh.x,us[:,i];
				title = string("surface at time t=",problem.times.ts[i]),
				label="")
	end
	if sav != []
		anim = @animate for i in range(1,length(ts))
			plotfig(i)
			#display(plt)
			#ylims!(plt[1,2],-1,18)
			#xlims!(plt[1,2],-10,75)
			#next!(prog)
		end
		gif(anim, string("anim-",sav,".gif"), fps=15); nothing
	end
	plotfig(length(ts))
	if sav != [] savefig(string("fig19f-",sav,".pdf")); end
	#
	# 	savefig(string("fig17-",sav,".pdf"));
	# 	p = plot(layout=(2,1))
	# 	fig_problem!(p,problem)
	# 	savefig(string("fig18-",sav,".pdf"));
	# 	save(problem,string("fig17-",sav));
	# end
end
figure19(0.9,false,sav="l09WGNdealias4")
