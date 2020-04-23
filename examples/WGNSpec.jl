export figspec

"""
Reproduces some figures in the work of C. Klein and V. Duchêne
"""
#using ShallowWaterModels
include("../src/dependencies.jl")


#---- Figure 14
"""
	figspec(c,SGN;sav)

Dispersive shock wave for SGN or WGN.
- Argument `c` is the velocity
- Solves SGN if `SGN` is `true`, WGN otherwise.
- Optional argument 'sav': if a string is given, then saves data and figures.

"""
function ell(a₀,a₁,k)
	h₀ = a₀
	h₂ = a₁ + a₀
	h₁ = (1-k^2)*a₁ + a₀
	return (h₀=h₀,h₁=h₁,h₂=h₂)
end

function figspecCW(p;SGN=true::Bool,sav=[],N=2^6,P=1)
	pa = ( μ  = 1, ϵ  = 1, N  = N, P=P )
	param=merge(p,pa)

	(η,u,v,parc,mesh)=CnoidalWaveWhithamGreenNaghdi(param,[0.];exact=true,SGN=true,P=P)
	@show parc

	param2=merge(pa,(L=P*parc.λ,))
	model = WhithamGreenNaghdi(param2;SGN=true)
	mesh2=Mesh(param2)
	plt = scatter()
	@showprogress 1 for ν=-π/parc.λ/P:π/parc.λ/P/1000:π/parc.λ/P
		L = Spectrum(model,η,u,parc.c,1im*ν)
		scatter!(plt,L.values,label="")
		# if maximum(real.(L.values))>0.001
		# 	return L,parc,η,u,v,mesh,ν
		# 	break
		# end
	end
	display(plt)
	return L,parc,η,u,v,mesh,mesh2

end

p=(h₀=0.4,h₁=1,h₂=1.2,ν=0)
p=merge(ell(0.3,0.1,0.99),(ν=1im*0.2-1im*π/λ/2,))

L,param,η,u,v,mesh,mesh2=figspecCW(ell(0.3,0.1,0.99);P=1,N=2^7)
L,param,η,u,v,mesh,ν=figspecCW(ell(0.3,0.1,0.999999999999999);P=1,N=2^8)
ylims!(0.2,0.25)
ylims!(0.2173064249,0.217306425)
scatter!(L.values)
λ=param.λ
plt = scatter()

@showprogress 1 for ν=-π/λ/2:π/λ/500:π/λ/2
	p=merge(ell(0.3,0.1,0.99),(ν=1im*ν,));
	L,=figspecCW(p;P=2,N=2^4);
	scatter!(plt,L.values,label="")
end
display(plt)
blabla

function figspecSW(c::Real,ν::Number;SGN=true::Bool,sav=[],N=2^5,L=8)
	param = ( μ  = 1, ϵ  = 1, δ = 0.01,
				N  = N,
	            L  = L,
				x0 = L/2,
				c=c,ν=ν)

	mesh = Mesh(param)
	model = WhithamGreenNaghdi(param;SGN=SGN, dealias = 0 )

	if c==0
		(η,u) = (0*mesh.x,0*mesh.x)
	elseif SGN==true
		function solη(x,c,ϵ,μ)
			(c^2-1)/ϵ*sech(sqrt(3*(c^2-1)/(c^2)/μ)/2*x)^2
		end
		η(c,x0)= solη.(mesh.x .-x0,c,param.ϵ,param.μ)
		u(c,x0)= c*η(c,x0)./(1 .+ param.ϵ*η(c,x0))
		(η,u) = (η(c,param.x0)+η(c,-param.x0),u(c,param.x0)-u(c,-param.x0))

	elseif SGN == false
		(η,u) = SolitaryWaveWhithamGreenNaghdi(
			mesh, param, η(c,param.x0);
			method=2, α = 1,
			tol =  1e-13, max_iter=10,
			iterative = false, q=1,
			verbose = true, SGN = false)
	end



	L = Spectrum(model,η,u,c,ν)

	plt = scatter(L.values,label="")
	display(plt)
	return L,param

end
ν=1im*0.1*π
a=(1-ν^2/3)^(1/2)
L,p=figspec(1.4,ν;N=2^9,L=32)
plt=scatter(L.values,label="")
scatter!([ν/a],label="a")
plot(Mesh(p).x,real.(L.vectors[1:Int(end/2),end-25].*exp.(1im*ν*Mesh(p).x)))
plt=plot()
#N=2^6;
#L,p=figspec(1.2,6/10*π/8;N=N,L=16)
#plt=scatter(L.values,label="0")
for ν = 0:π/8/10:π/8
	N=2^5
	L,p=figspec(10,ν;N=N,L=32)
	scatter!(plt,L.values[end-20:end],label=string(ν))
end
display(plt)
L,p=figspec(1.5,π/128/10;N=2^8,L=128)
test
function stab(L,param)
	model = WhithamGreenNaghdi(param;SGN=SGN, ktol=0, gtol=1e-12, iterate=true, precond = false, dealias = 0 )

	#
	# if sav != [] savefig(string("spectrum-c=",c,".pdf")); end
	# l = real.(L.values)
	# li = imag.(L.values)
	# l[abs.(li) .> 10] .= 0
	# k1=argmax(l)
	# l[k1]=0
	# k2=argmax(l)
	# k1 = 2*param.N-7
	# k2 = 2*param.N-7
	# η1 = L.vectors[1:Int(end/2),k1]
	# η2 = L.vectors[1:Int(end/2),k2]
	# δη = param.δ*real.(η1+η2)
	# v1 = L.vectors[Int(end/2)+1:end,k1]
	# v2 = L.vectors[Int(end/2)+1:end,k2]
	# δv = param.δ*real.(v1+v2)
	#
	# h = 1 .+ param.ϵ*η
	# DxF(v) = ifft(model.F₀ .* fft(v))
	# v= u - 1/3 ./h .* (DxF(h.^3 .*DxF(u)))
	#
	# init     = Init(mesh,η+δη,v+δv)
	# problem = Problem(model, init, param)
	#
	# solve!( problem )
	#
	# function fig(i)
	# 	ηinit = real.(ifft(problem.data.U[1][:,1]))
	# 	ηfin = real.(ifft(problem.data.U[i][:,1]))
	# 	plt = plot(layout=(1,2))
	# 	plot!(plt[1,1],fftshift(mesh.k),fftshift(log10.(abs.(fft(ηinit))));
	# 			title = "frequency",
	# 			label = "initial")
	# 	plot!(plt[1,1],fftshift(mesh.k),fftshift(log10.(abs.(fft(ηfin))));
	# 			label = "final")
	# 	plot!(plt[1,2],mesh.x,[ηinit ηfin];
	# 			title = string("error at time t=",problem.times.ts[i]),
	# 			label=["" ""])
	#
	# 	display(plt)
	# end
	#
	# ts = problem.times.ts
	#
	# (ηfin,vfin,ufin)   =  mapfro(model,last(problem.data.U))
	# (ηinit,vinit,uinit) = mapfro(model,first(problem.data.U))
	#
	# ηf(c,x0)= solη.(mesh.x .-x0,c,param.ϵ,param.μ)
	# ηinit = ηinit - ηf(c,0)
	# ηfin = ηfin - ηf(c,c*param.T)
	#
	# plt = plot(layout=(1,2))
	# plot!(plt[1,1],fftshift(mesh.k),fftshift(log10.(abs.(fft(ηinit))));
	# 		title = "frequency",
	# 		label = "initial")
	# plot!(plt[1,1],fftshift(mesh.k),fftshift(log10.(abs.(fft(ηfin))));
	# 		label = "final")
	# plot!(plt[1,2],mesh.x,[ηinit ηfin];
	# 		title = string("at time t=",ts[end]),
	# 		label=["zeta initial" "zeta final"])
	# # plot!(plt[1,2],mesh.x,[uinit ufin];
	# 		# title = string("at time t=",ts[end]),
	# 		# label=["u initial" "u final"])
	# display(plt)
	# if sav != [] savefig(string("final-c=",c,".pdf"));end
	#
	# x = mesh.x
	# k = mesh.k
	# us=zeros(param.N,length(ts));
	#
	# @showprogress 1 for i in 1:length(ts)
	# 	us[:,i].=real.(ifft(problem.data.U[i][:,1]))
	# end
	#
	# plt = plot()
	# plot!(plt,ts,maximum(abs.(us),dims=1)',
	# 		title="L infty norm",
	# 		label="")
	# display(plt)
	#
	# if sav != [] savefig(string("norm-c=",c,".pdf")); end
	#
	# # if sav != []
	# # 	anim = @animate for l in range(1,length(ts))
	# # 		fig(l)
	# # 	end
	# #
	# # 	gif(anim, string("anim-c=",c,".gif"), fps=15); nothing
	# # end
	#
	# fig(length(ts))
end
L=figspec(1.4,true)
scatter(L.values,label="")
