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
- `α=1.5` the parameter in the first initial data, describing the regularity,
- `(M,h₀=1/2,v₀=2) the parameters in the second initial data, describing the high frequency mode, the minimum depth, and the maximum velocity,
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
function IntegrateSV(;init=1,α=1.5,M=nothing,h₀=1/2,v₀=2,ϵ=1,L=π,N=2^9,T=0.5,dt =1e-5,dealias=true,smooth=false,Ns=nothing,label="Saint-Venant")
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
		if isnothing(M) M=(N÷3) end
		init = Init(x->(h₀-1)*cos.(x), x->2*sin.(x).+sin.(M*x)/M^2)
	else
		@error "argument init must be 1 or 2"
	end

	model = SaintVenant_fast(param; dealias = dealias, smooth = smooth, label = label) 
	solver=RK4(model.mapto(init))
	problem = Problem(model, init, param;solver=solver)
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

	problems=Problem[]
	E0=Float64[]
	E1=Float64[]
	
	if init == 1
		reference_problem=IntegrateSV(init=init,α=1.5,ϵ=1,L=π,N=2^15,T=0.05,dt = 1e-5,dealias=true,smooth=false,Ns=100,label="reference")
		Uref=reference_problem.data.U[end]/2^15

		for n=6:14
			p = IntegrateSV(init=init,α=1.5,ϵ=1,L=π,N=2^n,T=0.05,dt = 1e-5,dealias=true,smooth=smooth,Ns=1,label="N=2^$n")
			push!(problems,p)
			η,v,x=p.model.mapfro(p.data.U[end])
			U=[p.data.U[end][1] ; p.data.U[end][2]]/2^n
			Nmodes=length(x)÷2
			Ucomp=[ Uref[1][1:Nmodes];  Uref[1][end-Nmodes+1:end] ;
						Uref[2][1:Nmodes];  Uref[2][end-Nmodes+1:end]  ]


			push!(E0,norm(U-Ucomp,2)/norm(Ucomp,2))

			k=[sqrt.(1 .+ Mesh(x).k.^2);sqrt.(1 .+ Mesh(x).k.^2)]
			push!(E1,norm(k.*U-k.*Ucomp,2)/norm(k.*Ucomp,2))
		end
	elseif init == 2

		for n=14:-1:6
			p = IntegrateSV(init=init,h₀=0.5,v₀=2,ϵ=1,L=π,N=2^n,T=0.05,dt = 1e-5,dealias=true,smooth=smooth,Ns=1,label="N=2^$n")
			reference_problem=IntegrateSV(init=init,M=(2^n÷3),h₀=0.5,v₀=2,ϵ=1,L=π,N=2^15,T=0.05,dt = 1e-5,dealias=true,smooth=smooth,Ns=1,label="N=2^$n")

			Uref=reference_problem.data.U[end]/2^15
			η,v,x=p.model.mapfro(p.data.U[end])
			U=[p.data.U[end][1] ; p.data.U[end][2]]/2^n
			Nmodes=length(x)÷2
			Ucomp=[ Uref[1][1:Nmodes];  Uref[1][end-Nmodes+1:end] ;
						Uref[2][1:Nmodes];  Uref[2][end-Nmodes+1:end]  ]


			push!(E0,norm(U-Ucomp,2)/norm(Ucomp,2))

			k=[sqrt.(1 .+ Mesh(x).k.^2);sqrt.(1 .+ Mesh(x).k.^2)]
			push!(E1,norm(k.*U-k.*Ucomp,2)/norm(k.*Ucomp,2))
		end
	end
	return reference_problem,E0,E1

end


#--- Figures
plot_font = "Computer Modern"
default(fontfamily=plot_font)

#--- Experiment 1 : heap, sharp low-pass filter
reference_problem,E0,E1=Convergence(1;smooth=false)	
η0,v0,x=reference_problem.model.mapfro(reference_problem.data.U[1])
Fig1a=plot(x,[η0 v0],label=["\$\\eta^0\$" "\$u^0\$"],xlabel="\$x\$")
savefig(Fig1a,"Fig1a.pdf");savefig(Fig1a,"Fig1a.svg");

# anim = @animate for t in LinRange(0,reference_problem.times.tfin,101)
# 	plt = plot(reference_problem;var=[:surface,:fourier],T=t)
# 	ylims!(plt[1],(-0.5,1.1))
# end
# gif(anim, "anim1.gif", fps = 15)

Fig1b=plot(;xlabel="\$N\$",axis=:log)
Ns=2 .^(6:14);
scatter!(Fig1b,Ns,E1,label="\$E_1\$",color=1)
scatter!(Fig1b,Ns,E0,label="\$E_0\$",color=2)
plot!(Fig1b,Ns,Ns.^(-1),label="",color=1)
plot!(Fig1b,Ns,Ns.^(-2),label="",color=2)
savefig(Fig1b,"Fig1b.pdf");savefig(Fig1b,"Fig1b.svg");

#--- Experiment 2 : heap, smooth low-pass filter
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

#--- Experiment 3 : high Fourier mode, sharp low-pass filter
reference_problem,E0,E1=Convergence(2;smooth=false)	

η0,v0,x=reference_problem.model.mapfro(reference_problem.data.U[1])

Fig3a=plot(x,[η0 v0],label=["\$\\eta^0\$" "\$u^0\$"],xlabel="\$x\$")
plot!(x,v0.-2*sin.(x),
    xlim=(0,π),#ylim=(2-1e-2,2+1e-2),
    frame=:box,
    #yticks=false,
	color=:2,
    inset=bbox(0.6,0.6,0.35,0.25),
    subplot=2,
	label="\$u^0-2\\sin(x)\$"
)
savefig(Fig3a,"Fig3a.pdf");savefig(Fig3a,"Fig3a.svg");

# anim = @animate for t in LinRange(0,reference_problem.times.tfin,101)
# 	plt = plot(reference_problem;var=[:surface,:fourier],T=t)
# 	ylims!(plt[1],(-0.5,1.1))
# end
# gif(anim, "anim3.gif", fps = 15)

Fig3b=plot(;xlabel="\$N\$",axis=:log)
Ns=2 .^(6:14);
scatter!(Fig3b,Ns,E1[end:-1:1],label="\$E_1\$",color=1)
scatter!(Fig3b,Ns,E0[end:-1:1],label="\$E_0\$",color=2)
plot!(Fig3b,Ns,Ns.^(-1),label="",color=1)
plot!(Fig3b,Ns,Ns.^(-2),label="",color=2)
savefig(Fig3b,"Fig3b.pdf");savefig(Fig1b,"Fig3b.svg");


# Some experiments to play around experiment 3.
n=10
p = IntegrateSV(init=2,h₀=0.5,v₀=2,ϵ=1,L=π,N=2^n,T=0.05,dt = 1e-5,dealias=true,smooth=smooth,Ns=1,label="N=2^$n")
reference_problem=IntegrateSV(init=2,M=(2^n÷3),h₀=0.5,v₀=2,ϵ=1,L=π,N=2^12,T=0.05,dt = 1e-5,dealias=true,smooth=smooth,Ns=1,label="N=2^15")
plot(p,T=0,var=:fourier_velocity)
plot!(reference_problem,T=0,var=:fourier_velocity)
xlims!(-5,5)
η,v,x=solution(p;T=0)
ηref,vref,xref=solution(reference_problem;T=0)
plot(x,v.-2*sin.(x))
plot!(xref,vref.-2*sin.(xref))

plot(p;var=[:surface :velocity :fourier :fourier_velocity])
plot!(reference_problem;var=[:surface :velocity :fourier :fourier_velocity])

η,v,x=solution(p)
ηref,vref,xref=solution(reference_problem)

∂=1im*Mesh(x).k
diff(η;n=1)=real.(ifft(∂.^n.*fft(η)))

plot(x,diff(v;n=2),label="sharp")

∂=1im*Mesh(xref).k
diff(η;n=1)=real.(ifft(∂.^n.*fft(η)))
plot!(xref,diff(vref;n=2),label="sharp")


#--- Experiment 4 : instability in zero-depth situation
pb_smooth=IntegrateSV(init=2,h₀=0,v₀=2,ϵ=1,L=π,N=2^12,T=0.1,dt = 1e-5,dealias=true,smooth=true,Ns=1,label="smooth")
pb_sharp=IntegrateSV(init=2,h₀=0,v₀=2,ϵ=1,L=π,N=2^12,T=0.1,dt = 1e-5,dealias=true,smooth=smooth,Ns=1,label="sharp")
Fig4a=plot(pb_sharp,var=[:surface :velocity :fourier_velocity])
plot!(pb_smooth,var=[:surface :velocity :fourier_velocity])
savefig(Fig4a,"Fig4a.pdf");savefig(Fig4a,"Fig4a.svg");


η,v,x=solution(pb_sharp)
∂=1im*Mesh(x).k
diff(η;n=1)=real.(ifft(∂.^n.*fft(η)))

Fig4b=plot(x,diff(v;n=2),label="sharp")
ηs,vs,=solution(pb_smooth)
plot!(x,diff(vs;n=2),label="smooth")
title!("\$\\partial_x^2v\$")
xlabel!("\$x\$")
savefig(Fig4b,"Fig4b.pdf");savefig(Fig4b,"Fig4b.svg");

η0,v0,x=solution(pb_sharp;T=0)

norm(diff(v;n=2))./norm(diff(v0;n=2))
norm(diff(vs;n=2))./norm(diff(v0;n=2))
nothing
