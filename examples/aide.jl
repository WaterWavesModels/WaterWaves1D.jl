include("../src/dependencies.jl")

sav="fig10-long2"
problem = load(sav)
param=problem.param

function prepare(c::Real)

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

end


function fig(i)
	#(ηfin,vfin,ufin)   =  mapfro(model,problem.data.U[i])
	ηinit = real.(ifft(problem.data.U[1][:,1]))
	ηfin = real.(ifft(problem.data.U[i][:,1]))
	plt = plot(layout=(1,2))
	plot!(plt[1,1],fftshift(mesh.k),fftshift(log10.(abs.(fft(ηinit))));
			title = "frequency",
			label = "initial")
	plot!(plt[1,1],fftshift(mesh.k),fftshift(log10.(abs.(fft(ηfin))));
			label = "final")
	plot!(plt[1,2],mesh.x,[ηinit ηfin];
			title = string("error at time t=",problem.times.ts[i]),
			label=["" ""])
	ylims!(plt[1,2],-1,18)
	xlims!(plt[1,2],-10,75)

	display(plt)
end

ts = problem.times.ts
us=zeros(param.N,length(ts));
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

savefig(string("fig10-",sav,".pdf"));

anim = @animate for l in range(1,101)
	fig(l)
	#next!(prog)
end

gif(anim, "anim.gif", fps=15); nothing
