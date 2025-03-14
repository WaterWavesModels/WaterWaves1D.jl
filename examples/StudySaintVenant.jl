# #
# Reproduces the figures in
# the work of V. Duchêne and J. Marstrander
# on the numerical discretization of quasilinear systems
# #
export IntegrateSV,Figure
using WaterWaves1D,FFTW,Plots,LinearAlgebra,ProgressMeter;
#using JLD #(uncomment if using @save)

#--- Integration

"""
	IntegrateSV(;init,args)

Integrate in time the Saint-Venant (shallow water) system with an initial data depending on the provided `init`
- if `init=1`, then surface deformation `η(t=0,x)=1/2*exp(-|x|^α)*exp(-4x^2)` and velocity `v(t=0,x)=1/2*cos(x+1)*exp(-4x^2)` (with `α` provided as an optional argument, by default `α=3/2`)
- if `init=2`, then surface deformation `η(t=0,x)=exp(-x^2)` and velocity `v(t=0,x)=exp.(-x^2)*(sin(x)+sin(K*x)/K^2)` (with `K` provided as an optional argument, by default `K=100`)

Other arguments are optional:
- `α` the parameter in the first initial data,
- `ϵ` the nonlinearity parameter (default is `1`),
- `L` the half-length of the mesh (default is `π`),
- `N` the number of collocation points (default is `2^9`),
- `T` the final time of integration (default is `0.5`),
- `dt` the timestep (default is `1e-5`),
- `dealias`: no dealisasing if set to `0` or `false`, otherwise `1/(3*dealias)` modes are set to `0` (corresponding to standard 2/3 Orszag rule if `dealias` is set to `1` or `true`, which is default),
- `smooth`: A smooth low-pass filter (whose scaling is defined by ) if set to `0` or `false` (default), otherwise only `2/(3*dealias)*(1-smooth/2)` modes are kept untouched,
- `Ns` the number of stored computed times (default is all times),
- `label`: a label (string) for future use (default is ""Saint-Venant"").


Return `problem` of the solved problem 
"""
function IntegrateSV(;init=1,α=1.5,M=nothing,ϵ=1,L=π,N=2^9,T=0.5,dt =1e-5,dealias=true,smooth=false,Ns=nothing,label="Saint-Venant")
	if isnothing(Ns)
		param = ( ϵ  = ϵ,
				N  = N, L  = L,
	            T  = T, dt = dt )
	else
		param = ( ϵ  = ϵ,
				N  = N, L  = L,
	            T  = T, dt = dt,
				Ns = Ns )
	end

	mesh=Mesh(param)
	if init == 1
		init = Init(x->1/2*exp.(-abs.(x).^α).*exp.(-4*x.^2),x->1/2*cos.(x.+1).*exp.(-4*x.^2))
	elseif init == 2
		if isnothing(M) M=(N÷3)-1 end
		init = Init(x->1/2*cos.(2*x), x->sin.(2*x).+sin.(M*x)/M^2)
	else
		@error "argument init must be 1 or 2"
	end

	model = SaintVenant(param; dealias = dealias, smooth = smooth, label = label) 
	problem = Problem(model, init, param)
	solve!( problem )

	(ηfin,vfin)   =  solution(problem)

	# check energy preservation
	x=mesh.x;k=mesh.k;


	e(η,v) = 1/2* ( η.^2 .+ (1 .+ϵ*η).*v.^2);
	e(η1,v1,η2,v2) = 1/2* ( (η1-η2).*(η1+η2) .+ (1 .+ϵ*η2).*(v1-v2).*(v1+v2)
						.+ϵ*(η1-η2).*v1.^2 );
	η0=init.η(x);v0=init.v(x)
	error_energy = abs(sum(e(ηfin,vfin,η0,v0))/sum(e(η0,v0)))
	@info "normalized preservation of the total energy: $error_energy\n"

	return problem
end

#--- Experiments

"""
	Convergence(init;smooth)

Convergence rate associated with numerical experiments.
The Saint-Venant (shallow water) system is numerically integrated with several values for the number of collocation points.

- `init=1`: the initial data are chosen as for Figure 1 and 2
- `init=2`: the initial data are chosen as for Figure 3

If optional argument `smooth` is `true` (default is `false`): a smooth low-pass filter is used.

Return relevant plots and problems depending on the situation.
"""
function Convergence(init;smooth=false,name=nothing,anim=false)

	reference_problem=IntegrateSV(init=init,α=1.5,ϵ=1,L=π,N=2^15,T=0.05,dt = 1e-5,dealias=true,smooth=false,Ns=100,label="reference")
	Uref=reference_problem.data.U[end]/2^15

	problems=Problem[]
	E0=Float64[]
	E1=Float64[]
	for n=6:14
		p = IntegrateSV(init=init,α=1.5,ϵ=1,L=π,N=2^n,T=0.05,dt = 1e-5,dealias=true,smooth=smooth,Ns=1,label="N=2^$n")
		push!(problems,p)
		η,v,x=p.model.mapfro(p.data.U[end])
		U=p.data.U[end]/2^n
		Nmodes=length(x)÷2
		Ucomp=vcat(Uref[1:Nmodes,:],  Uref[end-Nmodes+1:end,:])

		push!(E0,norm(U-Ucomp,2)/norm(Ucomp,2))

		k=sqrt.(1 .+ Mesh(x).k.^2)
		push!(E1,norm(k.*U-k.*Ucomp,2)/norm(k.*Ucomp,2))
	end
	return reference_problem,E0,E1

end


#--- Figures
plot_font = "Computer Modern"
default(fontfamily=plot_font)

reference_problem,E0,E1=Convergence(1;smooth=false)	
η0,v0,x=reference_problem.model.mapfro(reference_problem.data.U[1])
Fig1a=plot(x,[η0 v0],label=["\$\\eta^0\$" "\$u^0\$"],xlabel="\$x\$")
savefig(Fig1a,"Fig1a.pdf");savefig(Fig1a,"Fig1a.svg");

anim = @animate for t in LinRange(0,reference_problem.times.tfin,101)
	plt = plot(reference_problem;var=[:surface,:fourier],T=t)
	ylims!(plt[1],(-0.5,1.1))
end
gif(anim, "anim.gif", fps = 15)

Fig1b=plot(;xlabel="\$N\$",axis=:log)
Ns=2 .^(6:14)
scatter!(Fig1b,Ns,E1,label="\$E_1\$",color=1)
scatter!(Fig1b,Ns,E0,label="\$E_0\$",color=2)
plot!(Fig1b,Ns,Ns.^(-1),label="",color=1)
plot!(Fig1b,Ns,Ns.^(-2),label="",color=2)
savefig(Fig1b,"Fig1b.pdf");savefig(Fig1b,"Fig1b.svg");


reference_problem,E0_smooth,E1_smooth=Convergence(1;smooth=true)	

Fig2=plot(;xlabel="\$N\$",axis=:log)
Ns=2 .^(6:14)
scatter!(Fig2,Ns,E1,label="\$E_1\$, sharp low-pass filter",color=1)
scatter!(Fig2,Ns,E0,label="\$E_0\$, sharp low-pass filter",color=2)
scatter!(Fig2,Ns,E1_smooth,label="\$E_1\$, smooth low-pass filter",color=3)
scatter!(Fig2,Ns,E0_smooth,label="\$E_0\$, smooth low-pass filter",color=4)

plot!(Fig2,Ns,Ns.^(-1),label="",color=1)
plot!(Fig2,Ns,Ns.^(-2),label="",color=2)

savefig(Fig2,"Fig2.pdf");savefig(Fig1b,"Fig2.svg");

# Table : experimental order of Convergence
round.(log2.(E0[1:end-1]./E0[2:end]),digits=2)
round.(log2.(E1[1:end-1]./E1[2:end]),digits=2)
round.(log2.(E0_smooth[1:end-1]./E0_smooth[2:end]),digits=2)
round.(log2.(E1_smooth[1:end-1]./E1_smooth[2:end]),digits=2)

nothing
