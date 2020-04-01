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

figure14(0.01,true;sav="SGN")

function sol(x,c,ϵ,μ)
	(c^2-1)/ϵ*sech.(sqrt(3*(c^2-1)/(c^2)/μ)/2*x).^2
end
ηGN(c)= sol(mesh.x,c,param.ϵ,param.μ)
uGN(c)= c*ηGN(c)./(1 .+ param.ϵ*ηGN(c))


#---- Figure 1
function figure1()
	param = ( μ  = 1,
			ϵ  = 1,
        	N  = 2^8,
            L  = 10*π,
						)
	mesh = Mesh(param)
	ηGN(c)= sol(mesh.x,c,param.ϵ,param.μ)
	uGN(c)= c*ηGN(c)./(1 .+ param.ϵ*ηGN(c))

	(η,u) = SolitaryWaveWhithamGreenNaghdi(
				mesh, merge(param,(c=1.1,)), ηGN(1.1);
				method=2, α = 0,
				tol =  1e-16, max_iter=4,
				ktol =1e-11, gtol = 1e-16,
				iterative = true, q=1,
				verbose = false, SGN = false)

	plt = plot(layout=(1,2))
	plot!(plt[1,1], mesh.x, [u uGN(1.1)];
	  title="c=1.1",
	  label=["WGN" "SGN"])

	plot!(plt[1,2], fftshift(mesh.k),
	  log10.(abs.(fftshift(fft(u))));
	  title="frequency",
	  label="WGN")
	savefig("fig1.pdf");
end

#---- Figure 2
function figure2()
	param = ( μ  = 1,
		ϵ  = 1,
      	N  = 2^10,
        L  = 10*π,)
	mesh = Mesh(param)
	ηGN(c)= sol(mesh.x,c,param.ϵ,param.μ)
	uGN(c)= c*ηGN(c)./(1 .+ param.ϵ*ηGN(c))

	(η,u) = SolitaryWaveWhithamGreenNaghdi(
				mesh, merge(param,(c=2,)), ηGN(2);
				method=2, α = 0,
				tol =  1e-12, max_iter=10,
				iterative = true, q=1,
				verbose = true, SGN = false)
	plt = plot(layout=(1,2))
	plot!(plt[1,1], mesh.x, [u uGN(2)];
	  title="c=2",
	  label=["WGN" "SGN"])
	plot!(plt[1,2], fftshift(mesh.k),
	  log10.(abs.(fftshift(fft(u))));
	  title="frequency",
	  label="WGN")
	savefig("fig2.pdf");

end
