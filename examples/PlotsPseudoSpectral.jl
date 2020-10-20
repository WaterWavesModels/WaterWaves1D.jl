# #
# Reproduces the figures in the work of B. Melinand and V. Duchêne
# on the Pseudo-spectral approximations of the water waves systems
# #
@info "Defines functions RWW2, StabilityExperiments etc."

using ShallowWaterModels,FFTW,Plots,LinearAlgebra,ProgressMeter;
include("../src/models/PseudoSpectral.jl")
include("../src/models/WaterWaves.jl")
include("../src/models/modifiedWW2.jl")
include("../src/models/ToyPseudoSpectral.jl")

include("../src/Figures.jl")
include("../src/LoadSave.jl")
using JLD

#----
"""
	`RWW2(;kwargs,name)`

Arguments are all optional:
- `eps` the nonlinearity parameter (default value `0.1`),
- `mu` the shallowness parameter (default value `1`);
- `delta` the regularisation parameter (default value `0.001`),
- `reg` the regularisation order (default value `2`),
- `dealias`: dealiasing if `1` (default), no dealiasing if `0`,
- `N` the number of modes (default value `2^16`),
- `T` final time of computation (default value `10`),
- `dt` timestep length (default value `0.01`),
- `L` half-length of the mesh (default value `10`),
- `pinit,vinit` initial data is `η₀(x)=exp(-|x|^pinit)`, `v₀(x)=vinit*x*exp(-|x|^pinit)` (default value `(2,1)`),
- `ns` data stored every `ns` computed step (default value `1`),
- `name`: a string used to save the figure as a name.pdf.

Return `(problem,plt)` where `problem` includes all computed data and `plt` the plot.
"""
function RWW2(;eps=0.1,mu=1,delta=0.001,reg=2,dealias=1,N=2^16,T=10,dt=0.01,L=20,pinit=2,vinit=0,ns=1,name=nothing)

	param = ( μ  = mu,
			ϵ  = eps,
        	N  = N, 	# number of collocation points
            L  = L,	# size of the mesh (-L,L)
            T  = T,		# final time of computation
            dt = dt, # timestep
			ns = ns,   	# data stored every ns computed step
				);
	@show((delta,reg,dealias))

	η₀(x) = exp.(-abs.(x).^pinit);
	v₀(x) = vinit*x.*exp.(-abs.(x).^pinit);

	init = Init(η₀,v₀);
	model=PseudoSpectral(param;order=2,dealias=dealias,δ=delta,reg=reg)
	problem=Problem(model, init, param)
	solve!( problem )

	plt = plot_solution(problem)
	display(plt)
	if name != nothing
		savefig(string(name,".pdf"));
		save(problem,name);
	end
	return problem,plt
end


#----
"""
	`StabilityExperiments(experiment=int,name=string)`

WW2 is numerically unstable, even with dealiasing. Regularizing restores stability.

- `int` is between 1 and 5
- `string` (optional) is used to save the figure as a string.pdf.

Return `(problem,plt)` where `problem` includes all computed data and `plt` the plot.
"""
function StabilityExperiments(;experiment,name=nothing)

	if experiment == 1
		@info "few modes without dealiasing (stable)"
		N = 2^9
		dealias = 0
		T = 10
		delta=0

	elseif experiment == 2
		@info "many modes without dealiasing (unstable)"
		N = 2^11
		dealias = 0
		T = 1.2
		delta=0

	elseif experiment == 3
		@info "few modes with dealiasing (stable)"
		N = 2^12
		dealias = 1
		T = 10
		delta = 0

	elseif experiment == 4
		@info "many modes with dealiasing (unstable)"
		N = 2^14
		dealias = 1
		T = 1.3
		delta = 0

	elseif experiment == 5
		@info "few modes with regularizing and without dealiasing (unstable)"
		N = 2^12
		dealias = 0
		T = 2.5
		delta = 0.01

	elseif experiment == 6
		@info "many modes with regularizing and dealiasing (stable)"
		N = 2^14
		dealias = 1
		T = 10
		delta = 0.01

	elseif experiment == 7
		@info "many modes with regularizing and dealiasing (stable)"
		N = 2^15
		dealias = 1
		T = 10
		delta = 0.002



	end

	pb,plt=RWW2(;eps=.1,mu=1,delta=delta,reg=2,dealias=dealias,N=N,T=T,dt=0.01,L=20,pinit=2,vinit=0,ns=1,name=name)
#	return pb,plt
end

"""
	`deltacritic(;delta_inf,delta_sup,tol_rel,kwargs)`

Seek the critical value for delta for existence of a solution.

Arguments are
- `delta_inf` initial lower value, `delta_sup` initial upper value
- `tol_rel` relative size of the interval
- `tol_abs` absolute size of the interval
- `kwargs` all arguments of `RWW2` (except `delta`)

Return `(delta_inf,delta_sup)` the computed interval.
"""
function deltacritic(;delta_inf=0,delta_sup=0.1,tol_rel=0.01,tol_abs=1e-10,eps=0.1,mu=1,reg=2,dealias=1,N=2^16,T=1,dt=0.01,L=20,pinit=2,vinit=0,ns=1,name=nothing)
	while (delta_sup-delta_inf) > max(tol_rel*delta_inf,tol_abs)
		delta = (delta_sup + delta_inf)/2
		@info string("Investigating the value of the critical delta on the interval ",@show (delta_inf,delta_sup))
		pb,plt=RWW2(;eps=eps,mu=mu,delta=delta,reg=reg,dealias=dealias,N=N,T=T,dt=dt,L=L,pinit=pinit,vinit=vinit,ns=ns,name=name)
		if any(isnan,pb.data.U[end][1])
			delta_inf = copy(delta)
		else
			delta_sup = copy(delta)
		end
	end
	@info delta_inf,delta_sup
	return delta_inf,delta_sup
end

"""
	`deltacheck(;delta_inf,delta_sup,kwargs)`

Checks whether the critical value for delta for existence of a solution is in a given interval.

Arguments are
- `delta_inf` initial lower value, `delta_sup` initial upper value
- `kwargs` all arguments of `RWW2` (except `delta`)

Return `(pb_inf,pb_sup)` the computed problem for the two values of delta.
"""
function deltacheck(;delta_inf=0,delta_sup=0.1,eps=0.1,mu=1,reg=2,dealias=1,N=2^16,T=1,dt=0.01,L=20,pinit=2,vinit=0)
	flag = []
	@info string("Checking whether the critical delta is on the interval ",@show (delta_inf,delta_sup))
	pb_inf,=RWW2(;eps=eps,mu=mu,delta=delta_inf,reg=reg,dealias=dealias,N=N,T=T,dt=dt,L=L,pinit=pinit,vinit=vinit)
	if !any(isnan,pb_inf.data.U[end][1])
		push!(flag,"delta_inf is too large")
	end
	pb_sup,=RWW2(;eps=eps,mu=mu,delta=delta_sup,reg=reg,dealias=dealias,N=N,T=T,dt=dt,L=L,pinit=pinit,vinit=vinit)
	if any(isnan,pb_sup.data.U[end][1])
		push!(flag,"delta_sup is too small")
	end
	if flag == []
		@info "success!"
	else
		@error flag
	end
	return flag,pb_inf,pb_sup
end

# deltacheck(eps=0.1,delta_inf=0.000998,delta_sup=0.001008,T=1,N=2^16,vinit=1)
# deltacheck(eps=0.1,delta_inf=0.00148,delta_sup=0.00152,T=5,N=2^16,vinit=1)
# deltacheck(eps=0.1,delta_inf=0.00148,delta_sup=0.00152,T=10,N=2^16,vinit=1)
# deltacheck(eps=0.1,delta_inf=0.00148,delta_sup=0.00152,T=15,N=2^16,vinit=1)

"""
	`CriticalDelta(Eps;kwargs)`

Plot critical values for delta corresponding to several values of epsilon

Arguments are
- `Eps` values of epsilon
- `kwargs` all arguments of `RWW2` (except `delta`)

Return `(Eps,Delta_inf,Delta_sup,plt)`
"""
function CriticalDelta(Eps;T=1,tol_rel=0.01,tol_abs=1e-10,mu=1,reg=2,dealias=1,N=2^16,dt=0.01,L=20,pinit=2,vinit=0,name=nothing)
	Delta_inf = empty(Eps); Delta_sup = empty(Eps); Eps_comp = empty(Eps);
	param = (tol_abs=tol_abs,tol_rel=tol_rel,mu=mu,reg=reg,dealias=dealias,N=N,T=T,dt=dt,L=L,pinit=pinit,vinit=vinit)
	for eps in Eps
		(delta_inf,delta_sup)=deltacritic(eps=eps,delta_inf=0,delta_sup=eps^2,tol_rel=tol_rel,tol_abs=tol_abs,mu=mu,reg=reg,dealias=dealias,N=N,T=T,dt=dt,L=L,pinit=pinit,vinit=vinit)
		push!(Delta_inf,delta_inf);push!(Delta_sup,delta_sup);push!(Eps_comp,eps);
		if name != nothing
			save(string(name,".jld"),"Eps_comp",Eps,"Delta_inf",Delta_inf,"Delta_sup",Delta_sup,"param",param);
		end
	end
	plt=scatter(log10.(Eps_comp),log10.(Delta_inf))
	scatter!(plt,log10.(Eps_comp),log10.(Delta_sup))
	xlabel!(plt,"log(epsilon)")
	ylabel!(plt,"log(delta)")
	display(plt)
	if name != nothing
		savefig(string(name,".pdf"));
	end

	return Eps_comp,Delta_inf,Delta_sup,plt
end

"""
	`CheckCriticalDelta(Eps,Delta_inf,Delta_sup;kwargs)`

Checks critical value intervals for several values of epsilon

Arguments are
- `Eps` values of epsilon,
- `Delta_inf` lower values for delta,
- `Delta_sup` upper values
- `kwargs` all arguments of `RWW2` (except `delta`)

Return `flag`
"""
function CheckCriticalDelta(Eps,Delta_inf,Delta_sup;T=1,mu=1,reg=2,dealias=1,N=2^16,dt=0.01,L=20,pinit=2,vinit=0,name=nothing)
	param = (mu=mu,reg=reg,dealias=dealias,N=N,T=T,dt=dt,L=L,pinit=pinit,vinit=vinit)
	Delta_inf_comp = empty(Eps); Delta_sup_comp = empty(Eps); Eps_comp = empty(Eps); flag_comp=[];
	for i in 1:length(Eps)
		flag,=deltacheck(eps=Eps[i],delta_inf=Delta_inf[i],delta_sup=Delta_sup[i],mu=mu,reg=reg,dealias=dealias,N=N,T=T,dt=dt,L=L,pinit=pinit,vinit=vinit)
		push!(Delta_inf_comp,Delta_inf[i]);push!(Delta_sup_comp,Delta_sup[i]);push!(Eps_comp,Eps[i]);
		if flag == []
			push!(flag_comp,"OK");
		else
			push!(flag_comp,flag);
		end
		@show flag_comp
		if name != nothing
			save(string(name,".jld"),"flag",flag_comp,"Eps",Eps_comp,"Delta_inf",Delta_inf_comp,"Delta_sup",Delta_sup_comp,"param",param);
		end
	end

	return flag_comp
end

# E,Di,Ds,plt=PlotCriticalDelta(0.03:0.01:0.1;vinit=1,name="testpct")

"""
	`PlotPrecision(;delta_inf,delta_sup,delta_step,kwargs)`

Compares solutions for several values of delta

Arguments are
- `delta_inf,delta_sup,delta_step`: deltas are in `Delta=range(delta_inf,stop=delta_sup,step=delta_step)`
- `kwargs` all arguments of `RWW2` (except `delta`)

Return `(Delta,Diff)` where `Diff` is the l^infty norm of the difference between the solution with and without regularization
"""
function PlotPrecision(;delta_inf=0.01,delta_sup=0.1,delta_length=10,T=1,eps=0.1,mu=1,reg=2,dealias=1,N=2^12,dt=0.01,L=20,pinit=2,vinit=0,name=nothing)

	param = ( μ  = mu,
			ϵ  = eps,
			N  = N, 	# number of collocation points
			L  = L,	# size of the mesh (-L,L)
			T  = T,		# final time of computation
			dt = dt, # timestep
			ns = 1,   	# data stored every ns computed step
				);

	η₀(x) = exp.(-abs.(x).^pinit);
	v₀(x) = vinit*x.*exp.(-abs.(x).^pinit);

	init = Init(η₀,v₀);
	model=WaterWaves(param)
	model1=modifiedWW2(param;dealias=dealias)
	model0=PseudoSpectral(param;order=2,δ=0,dealias=dealias,reg=reg)
	#model3=PseudoSpectral(param;order=3,δ=0,dealias=dealias,reg=reg)
	#model4=PseudoSpectral(param;order=2,δ=0.01,dealias=dealias,reg=reg)

	pb=Problem(model, init, param)
	pb1=Problem(model1, init, param)
	pb0=Problem(model0, init, param)
	#pb3=Problem(model3, init, param)
	#pb4=Problem(model4, init, param)


	solve!( pb )
	solve!( pb1 )
	solve!( pb0 )
	#solve!( pb3 )
	#solve!( pb4 )

	(η,v,x)=solution(pb)
	(η1,v1)=solution(pb1,x=x);E1=norm(η1-η,Inf);
	(η0,v0)=solution(pb0,x=x);E0=norm(η0-η,Inf);


	Delta = 10 .^(range(log10.(delta_inf),stop=log10.(delta_sup),length=delta_length))
	Diff = similar(Delta)
	for i in 1:length(Delta)
		pbdelta,=RWW2(;eps=eps,mu=mu,delta=Delta[i],reg=reg,dealias=dealias,N=N,T=T,dt=dt,L=L,pinit=pinit,vinit=vinit)
		ηd,vd = solution(pbdelta,x=x)
		Diff[i] = norm(ηd-η,Inf)
		if name != nothing
			save(string(name,".jld"),"Delta",Delta,"Diff",Diff,"param",param);
		end
	end
	plt=scatter(Delta,Diff,yscale=:log10,xscale=:log10,label="")
	xlabel!(plt,"\$\\delta\$ (in log scale)")
	ylabel!(plt,"error (in log scale)")
	display(plt)
	return Delta,Diff,E0,E1,plt
end

function PlotFig4(p;name=nothing)
	param=(pinit=p,delta_inf=0.01,delta_sup=1,delta_length=20,N=2^12,dt=0.01,T=10,mu=1,reg=2,dealias=1,L=20,epsa=0.1,epsb=0.05,epsc=0.2)
	Delta,Diffa = PlotPrecision(pinit=p,delta_inf=0.01,delta_sup=1,delta_length=20,N=2^12,dt=0.01,T=10,eps=0.1,mu=1,reg=2,dealias=1,L=20);
	Delta,Diffb = PlotPrecision(pinit=p,delta_inf=0.01,delta_sup=1,delta_length=20,N=2^12,dt=0.01,T=10,eps=0.05,mu=1,reg=2,dealias=1,L=20);
	Delta,Diffc = PlotPrecision(pinit=p,delta_inf=0.01,delta_sup=1,delta_length=20,N=2^12,dt=0.01,T=10,eps=0.2,mu=1,reg=2,dealias=1,L=20);

	plt=plot()
	scatter!(plt,Delta,Diffa,label="\$\\epsilon=0.01\$",legend=:topleft)
	scatter!(plt,Delta,Diffb,label="\$\\epsilon=0.05\$",legend=:topleft,marker=:d)
	scatter!(plt,Delta,Diffc,label="\$\\epsilon=0.2\$",legend=:topleft,marker=:r)
	if name !=nothing
		savefig(plt,string(name,".pdf"))
		savefig(plt,string(name,"-log.svg"))
	end

	plt=plot()
	scatter!(plt,Delta,Diffa,yscale=:log10,xscale=:log10,label="\$\\epsilon=0.01\$",legend=:topleft)
	scatter!(plt,Delta,Diffb,yscale=:log10,xscale=:log10,label="\$\\epsilon=0.05\$",legend=:topleft,marker=:d)
	scatter!(plt,Delta,Diffc,yscale=:log10,xscale=:log10,label="\$\\epsilon=0.2\$",legend=:topleft,marker=:r)
	if name !=nothing
		savefig(plt,string(name,"-log.pdf"))
		savefig(plt,string(name,"-log.svg"))

		save(string(name,".jld"),"Delta",Delta,"Diffa",Diffa,"Diffb",Diffb,"Diffc",Diffc,"param",param);

	end

	display(plt)
	return Delta,Diffa,Diffb,Diffc
end

#PlotFig4(1,name="precision1")
#PlotFig4(3,name="precision3")
#PlotFig4(2,name="precision2")

function toadd()
	η,v = solution(pb0)


	param = (eps=eps,mu=mu,reg=reg,dealias=dealias,N=N,T=T,dt=dt,L=L,pinit=pinit,vinit=vinit)

	Delta = range(delta_inf,stop=delta_sup,step=delta_step)
	Diff = similar(Delta)
	for i in 1:length(Delta)
		pbdelta,=RWW2(;eps=eps,mu=mu,delta=Delta[i],reg=reg,dealias=dealias,N=N,T=T,dt=dt,L=L,pinit=pinit,vinit=vinit)
		ηd,vd = solution(pbdelta)
		Diff[i] = norm(η-ηd,Inf)
		if name != nothing
			save(string(name,".jld"),"Delta",Delta,"Diff",Diff,"param",param);
		end
	end
	plt=scatter(log10.(Delta),log10.(Diff),label="")
	xlabel!(plt,"log(delta)")
	ylabel!(plt,"log(error)")
	display(plt)
	return Delta,Diff,plt
end




# #### Some experiments for PlotPrecision
# Delta2,Diff2=PlotPrecision(eps=0.1,delta_inf=0,delta_sup=0.1,delta_step=0.001,T=10,N=2^12,dt=0.01)
# #Delta2,Diff2=PlotPrecision(eps=0.1,delta_inf=0,delta_sup=0.1,delta_step=0.001,T=10,N=2^12,dt=0.001)
# Delta1,Diff1=PlotPrecision(eps=0.1,delta_inf=0,delta_sup=0.1,delta_step=0.001,T=10,N=2^11,dt=0.01,pinit=1)
# Delta3,Diff3=PlotPrecision(eps=0.1,delta_inf=0,delta_sup=0.1,delta_step=0.001,T=10,N=2^11,dt=0.01,pinit=3)
#
# scatter(Delta1,log10.(Diff3),label="")
# scatter(log10.(Delta2),log10.(Diff2),label="")
# scatter!(log10.(Delta2[Delta2.<0.0015]),log10.(Diff2[Delta2.<0.0015]),label="")
# plot!(log10.(Delta3),3*log10.(Delta3),label="slope 3",legend=:topleft)
# plot!(log10.(Delta3),0*log10.(Delta3).+log10(0.0009998292076658089),label="typical numerical error",legend=:topleft)
# xlabel!("log(delta)")
# ylabel!("log(error)")
#
#
# deltacheck(;eps=0.1,delta_inf=0.004825,delta_sup=0.00485,T=5,N=2^16,pinit=1)
# deltacheck(eps=0.1,delta_inf=0.0015,delta_sup=0.00151,T=5,N=2^16,pinit=2);
# deltacheck(;eps=0.1,delta_inf=0.00215,delta_sup=0.002175,T=5,N=2^16,pinit=3)
#
# pb1,=RWW2(delta=0.01,eps=0.1,N=2^11,T=1,dt=0.01,pinit=1)
# pb2,=RWW2(delta=0.01,eps=0.1,N=2^16,T=1,dt=0.01,pinit=1)
#
# η1,v1,x1 = solution(pb1)
# η2,v2,x2 = solution(pb2;x=x1)
# @show norm(η1-η2,Inf)
# # pinit= 3 Error:6.778345691960119e-7,8.547625869881337e-7,8.980761364574263e-7
# 0.0009998292076658089
# #### end


##### Pour voir croissance en temps des haut modes
# ts=pb.times.ts;		Ns=pb.times.Ns
# Nη=[];Nv=[];
# for t in ts
# 	(η,v)=solution(pb,t=t)
# 	push!(Nη,maximum(log10.(abs.(fft(η)[250:end-250]))))
# 	push!(Nv,maximum(log10.(abs.(fft(v)[250:end-250]))))
# end
# plot(ts,[Nη Nv]);plot!(ts,16*ts.-15)
# ylims!(-15,0)
#
#
# A=[]
# B=[]
# for t=0:0.01:2
# 	(η,v)=solution(pb,t=t)
# 	push!(B,maximum(log.(abs.(fftshift(fft(η)[200:end-200])))))
# 	push!(A,t)
# end
# plot(A,30*A.-48);plot!(A,B)
#
# t=1
# cutoff = k -> (1-exp(-1/abs(k)^2))^(reg/2)
# δ=pb.model.kwargs.δ
# Π = cutoff.(δ*k)
# k=fftshift(pb.mesh.k)
# (η,v)=solution(pb,t=t)
# plot(k,[log.(abs.(fftshift(fft(η)))).+48 sqrt.(abs.(k).^(3/2).*Π.^2).^(2/3)*1.2*t])
#####
nothing


function TestBonaetal()
	param = ( μ  = Inf,
		ϵ  = 1,
		N  = 320*8, 	# number of collocation points
		L  = π,	# size of the mesh (-L,L)
		T  = 10,		# final time of computation
		dt = 0.001, # timestep
		ns = 1,   	# data stored every ns computed step
			);

	J=32

	η₀(x) = zero(x);
	v₀(x) = sin.(J*x)/J^(1/2)/8;

	init = Init(η₀,v₀);
	model=PseudoSpectral(param;order=2,dealias=1,δ=0*J/320/2/2/2,reg=2)
	problem=Problem(model, init, param)
	solve!( problem )

	plt = plot_solution(problem,t=10,velocity=true)
end

function WhatHappensAtCriticalDelta()
	pb,plt=StabilityExperiments(experiment=7)

	function Norm(;t=1)
		(η,v)=solution(pb,t=t)
		Dv=real.(ifft(fft(v)./(1im*k.+eps(1.)).*k.*tanh.(k)))
		plot(x,[v Dv])
		return sum(Dv.*Dv)
	end
	N=empty([.1])
	for t in pb.times.ts
		push!(N,Norm(;t=t))
	end

	function Norm(;i=1)
		fftη=pb.data.U[i][:,1][250:end-250]
		#plot(k[250:end-250],abs.(fftη).+eps(1.),yscale=:log10)
		return norm(fftη,2)
	end
	N2=empty([.1])
	@progress 1 for i in 1:length(pb.times.ts)
		push!(N2,Norm(;i=i))
	end

	plot(pb.times.ts,[ normalize(N) normalize(N2)])
end


function PlotFig2()
	Delta_inf = [2.8125e-5, 0.0001, 0.00018984375000000002, 3.9374999999999995e-5, 5.9765625e-5, 0.00014160156250000002, 0.0003062500000000001, 0.00042968750000000006, 0.0006436523437500001, 0.0009, 0.0011718750000000002, 0.0015187500000000001, 0.0020214843750000003, 0.002890625000000001, 0.003662109375, 0.0045703125, 0.006250000000000001, 0.008056640625, 0.010727539062500001, 0.015000000000000003, 0.0380859375]
	Delta_sup = [2.890625e-5, 0.000103125, 0.000196875, 4.0499999999999995e-5, 6.15234375e-5, 0.00014648437500000002, 0.00031875, 0.0004492187500000001, 0.0006601562500000002, 0.000925, 0.0012109375000000002, 0.001575, 0.002109375, 0.002968750000000001, 0.0037841796875, 0.0047460937499999994, 0.006562500000000001, 0.00830078125, 0.011140136718750002, 0.015625, 0.0390625]
	Eps = [0.01, 0.02, 0.03, 0.012, 0.015, 0.025, 0.04, 0.05, 0.065, 0.08, 0.1, 0.12, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.65, 0.8, 1.0]
	Delta2_sup = [2.890625e-5, 4.0499999999999995e-5, 6.15234375e-5, 0.000103125, 0.00014648437500000002, 0.000196875, 0.00031250000000000006, 0.0004492187500000001, 0.0006271484375000001, 0.0008500000000000001, 0.0011328125000000001, 0.0014625, 0.00193359375, 0.0025781250000000006, 0.0032958984375, 0.0040429687500000006, 0.005312500000000001, 0.006591796875, 0.008045654296875, 0.009687500000000002, 0.0244140625]
	Delta2_inf = [2.8125e-5, 3.9374999999999995e-5, 5.9765625e-5, 0.0001, 0.00014160156250000002, 0.00018984375000000002, 0.00030000000000000003, 0.00042968750000000006, 0.0006106445312500001, 0.000825, 0.00109375, 0.00140625, 0.0018457031250000001, 0.0025000000000000005, 0.003173828125, 0.0038671875, 0.005156250000000001, 0.00634765625, 0.00783935546875, 0.009375000000000001, 0.0234375]
	Eps2 = [0.01, 0.012, 0.015, 0.02, 0.025, 0.03, 0.04, 0.05, 0.065, 0.08, 0.1, 0.12, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.65, 0.8, 1.0]
	scatter(Eps,Delta_sup/2+Delta_inf/2,yerr=Delta_sup/2-Delta_inf/2,marker =  stroke(3),markersize=0,yscale=:log10,xscale=:log10,legend=:none,label="")
	scatter!(Eps2,Delta2_sup/2+Delta2_inf/2,yerr=Delta2_sup/2-Delta2_inf/2,marker =  stroke(3),markersize=0,yscale=:log10,xscale=:log10,legend=:none,label="")
	plot!(Eps,Eps.^2,label="\$\\delta=\\epsilon^2\$",legend=:topleft)
	xlabel!("\$\\epsilon\$")
	ylabel!("\$\\delta\$")
end

function blowup(pb)
	T = try pb.times.ts[indexin(NaN,real.(sum.(pb.data.U)))][1]
	catch
		@warn "no blowup"
		T=last(pb.times.ts)
	end
end
#
# global T=[]
# for eps in [0.8 0.4 0.2 0.1 0.05 0.025]
# 	pb,=RWW2(T=5,eps=eps)
# 	push!(T,blowup(pb))
# 	@show T
# end
#
# T= [0.45, 0.65, 1.0, 1.9, 5.0, 5.0]

#pb,=RWW2(T=5,eps=0.1,delta=0.01,vinit=-1)

pb1d2,= RWW2(T=10,reg=1,delta=0.01,pinit=3,N=2^16)
pb34d2,= RWW2(T=10,reg=3/4,delta=0.01,pinit=3,N=2^16)
pb12d2,= RWW2(T=10,reg=1/2,delta=0.01,pinit=3,N=2^16)
pb14d2,= RWW2(T=10,reg=1/4,delta=0.01,pinit=3,N=2^16)

pb1d3,= RWW2(T=10,reg=1,delta=0.01,pinit=3,N=2^18)
pb34d3,= RWW2(T=10,reg=3/4,delta=0.01,pinit=3,N=2^18)
pb12d3,= RWW2(T=10,reg=1/2,delta=0.01,pinit=3,N=2^18)
pb14d3,= RWW2(T=10,reg=1/4,delta=0.01,pinit=3,N=2^18)


pb1d1,= RWW2(T=10,reg=1,delta=0.1,pinit=3,N=2^16)
pb34d1,= RWW2(T=10,reg=3/4,delta=0.1,pinit=3,N=2^16)
pb12d1,= RWW2(T=10,reg=1/2,delta=0.1,pinit=3,N=2^16)
pb14d1,= RWW2(T=10,reg=1/4,delta=0.1,pinit=3,N=2^16)


pb1d0,= RWW2(T=10,reg=1,delta=0.1,pinit=3,N=2^18)
pb34d0,= RWW2(T=10,reg=3/4,delta=0.1,pinit=3,N=2^18)
pb12d0,= RWW2(T=10,reg=1/2,delta=0.1,pinit=3,N=2^18)
pb14d0,= RWW2(T=10,reg=1/4,delta=0.1,pinit=3,N=2^18)

global T=[]
for reg=1/4-0.05:0.01:1/2
	pb, =RWW2(T=5,reg=reg,delta=0.1,pinit=3,N=2^18)
	push!(T,blowup(pb))
end

# delta=0.1,pinit=3,N=2^18
Reg = 1/4:0.01:1/2
T= [0.79, 0.82, 0.85, 0.88, 0.92, 0.96, 1.0, 1.05, 1.09, 1.14, 1.2, 1.24, 1.3, 1.35, 1.41, 1.48, 1.54, 1.62, 1.7, 1.8, 1.91, 2.0, 2.0, 2.0, 2.0, 2.0]
plot(Reg,T)

# delta=0.01,pinit=3,N=2^18
Reg = 1/4:0.01:1/2
T= [0.79, 0.82, 0.85, 0.88, 0.92, 0.96, 1.0, 1.05, 1.09, 1.14, 1.2, 1.24, 1.3, 1.35, 1.41, 1.48, 1.54, 1.62, 1.7, 1.8, 1.91, 2.0, 2.0, 2.0, 2.0, 2.0]
plot(Reg,T)


pb1, =RWW2(T=20,reg=2,eps=0.1,delta=0.1,pinit=3,N=2^10,dt=0.01)
pb32, =RWW2(T=20,reg=2,eps=0.1,delta=0.03,pinit=3,N=2^10,dt=0.01)
pb2, =RWW2(T=20,reg=2,eps=0.1,delta=0.01,pinit=3,N=2^10,dt=0.01)


plot(pb1.times.ts,norm.(pb32.data.U-pb1.data.U))
plot(pb1.times.ts,norm.(pb2.data.U-pb32.data.U))
plot(pb1.times.ts,norm.(pb2.data.U-pb1.data.U))

global N=[]
function Norm(pb)
	for u in pb32.data.U
		push!(N,norm(u .*k.^30))
	end
end
#
# plot_solution(pb1,label="")
# savefig("order1.pdf")
# savefig("order1.svg")
# plot_solution(pb34,label="")
# savefig("order34.pdf")
# savefig("order34.svg")
# plot_solution(pb12,label="")
# savefig("order12.pdf")
# savefig("order12.svg")
# plot_solution(pb14,label="")
# savefig("order14.pdf")
# savefig("order14.svg")

function RWW2b(;J=100,s=3,eps=0.1,mu=1,delta=0.001,reg=2,dealias=1,N=2^16,T=10,dt=0.01,L=20,pinit=2,vinit=0,ns=1,name=nothing)

	param = ( μ  = mu,
			ϵ  = eps,
        	N  = N, 	# number of collocation points
            L  = L,	# size of the mesh (-L,L)
            T  = T,		# final time of computation
            dt = dt, # timestep
			ns = ns,   	# data stored every ns computed step
				);
	@show((delta,reg,dealias))

	η₀(x) = exp.(-abs.(x).^pinit).*(1 .+ cos.(J*x)/J^s);
	v₀(x) = exp.(-abs.(x).^pinit);

	init = Init(η₀,v₀);
	model=PseudoSpectral(param;order=2,dealias=dealias,δ=delta,reg=reg)
	problem=Problem(model, init, param)
	solve!( problem )

	plt = plot_solution(problem)
	display(plt)
	if name != nothing
		savefig(string(name,".pdf"));
		save(problem,name);
	end
	return problem,plt
end

pb1, =RWW2b(J=700,s=3,T=2,reg=2/8,eps=0.1,delta=0.01,pinit=2,vinit=1,N=2^14,dt=0.01)

function Norm(pb,i)
	maximum(abs.(pb.data.U[i][end÷4:3*end÷4,1]))
end

N=[]
for t in 1:pb1.times.Ns
	push!(N,Norm(pb1,t))
end
