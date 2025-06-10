# #
# Reproduces the figures in
# the work of V. Duchêne and J. Marstrander
# on the numerical discretization of quasilinear systems
# #
export IntegrateSV,IntegrateSV2D,Convergence,Convergence2D,Figure1,Figure2,Figure3,Figure4,Figure5,Figure6,Figure7
using WaterWaves1D,FFTW,Plots,LinearAlgebra,ProgressMeter;

plot_font = "Computer Modern"
default(fontfamily=plot_font)


#------------------- 
#--- 1D PROBLEMS ---
#-------------------


#--- Integration

"""
	IntegrateSV(;init,args)

Integrate in time the Saint-Venant (shallow water) system with an initial data depending on the provided `init`
- if `init=1`, then surface deformation `η(t=0,x)=1/2*exp(-|x|^α)*exp(-4x^2)` and velocity `v(t=0,x)=0` (with `α` provided as an optional argument, by default `α=3/2`)
- if `init=2`, then surface deformation `η(t=0,x)=exp(-x^2)` and velocity `v(t=0,x)=exp.(-x^2)*(sin(x)+sin(N*x)/N^2)` (with `N` provided as an optional argument)

Other arguments are optional:
- `α` the parameter in the first initial data, prescribing the regularity (default is `3/2`),
- `(N,h₀,v₀) the parameters in the second initial data, describing the high frequency mode, the minimum depth, and the maximum velocity 
(by default, `h₀=1/2,v₀=2` and `N` is 2/3 of the number of collocation points),
- `ϵ` the nonlinearity parameter (default is `1`),
- `L` the half-length of the mesh (default is `π`),
- `M` the number of collocation points (default is `2^9`),
- `T` the final time of integration (default is `0.5`),
- `dt` the timestep (default is `1e-5`),
- `dealias`: no dealisasing if set to `0` or `false`, otherwise `1/(3*dealias)` modes are set to `0` (corresponding to standard 2/3 Orszag rule if `dealias` is set to `1` or `true`, which is default),
- `smooth`: A smooth low-pass filter (whose scaling is defined by `N`) if set to `0` or `false` (default), otherwise only `2/(3*dealias)*(1-smooth/2)` modes are kept untouched,
- `Ns` the number of stored computed times (default is all times),
- `label`: a label (string) for future use (default is ""Saint-Venant"").


Return `problem` of the solved problem 
"""
function IntegrateSV(;init=1,α=1.5,N=nothing,h₀=1/2,v₀=2,ϵ=1,L=π,M=2^9,T=0.5,dt =1e-5,dealias=true,smooth=false,Ns=nothing,label="Saint-Venant")
	if isnothing(Ns)
		param = ( ϵ  = ϵ,
				N  = M, L  = L,
	            T  = T, dt = dt )
	else
		param = ( ϵ  = ϵ,
				N  = M, L  = L,
	            T  = T, dt = dt,
				Ns = Ns )
	end

	mesh=Mesh(param)
	if init == 1
		init = Init(x->1/2*exp.(-abs.(x).^α).*exp.(-4*x.^2),x->1/2*cos.(x.+1).*exp.(-4*x.^2))
	elseif init == 2
<<<<<<< Updated upstream
		if isnothing(M) M=(N÷3) end
		init = Init(x->(h₀-1)*cos.(x), x->2*sin.(x).+sin.(M*x)/M^2)
	elseif init == 3
		if isnothing(M) M=(N÷3) end
		init = Init(x->(h₀-1)*cos.(x).-sin.(M*x)/M^2 .+cos.((M-1)*x)/M^2, x->2*sin.(x).+sin.(M*x)/M^2 .-cos.((M-1)*x)/M^2)
=======
		if isnothing(N) N=(M÷3) end
		init = Init(x->(h₀-1)*cos.(x), x->v₀ .* sin.(x).+sin.(N*x)/N^2)
	elseif init == 3
		f(x)=1/2*exp.(-4*x.^2)
		init = Init(x->f(x),x->2*sqrt.(1 .+ϵ*f(x)).-2)
>>>>>>> Stashed changes
	else
		@error "argument init must be 1 or 2"
	end

	model = SaintVenant_fast(param; dealias = dealias, smooth = smooth, label = label) 
	solver=RK4(model.mapto(init))
	problem = Problem(model, init, param;solver=solver)
	solve!( problem )


	# check energy preservation
	e(η,v) = 1/2* ( η.^2 .+ (1 .+ϵ*η).*v.^2);
	e(η1,v1,η2,v2) = 1/2* ( (η1-η2).*(η1+η2) .+ (1 .+ϵ*η2).*(v1-v2).*(v1+v2)
						.+ϵ*(η1-η2).*v1.^2 );
	
	(ηfin,vfin)   =  solution(problem)
	(η0,v0)=(init.η(mesh.x),init.v(mesh.x))

	error_energy = abs(sum(e(ηfin,vfin,η0,v0))/sum(e(η0,v0)))
	@info "normalized preservation of the total energy: $error_energy\n"

	return problem
end


"""
	Convergence(init;smooth,T=0.5)

Convergence rate associated with numerical experiments.
The Saint-Venant (shallow water) system is numerically integrated with several values for the number of collocation points.

- `init=1`: the initial data are chosen as for Figure 1 and 2 (heap of water)
- `init=2`: the initial data are chosen as for Figure 3 (high-frequency component)

If optional argument `smooth` is `true` (default is `false`): a smooth low-pass filter is used.
The optional argument `T` is the final time of computation (default is `0.5`).

Return the reference problem and relative errors.
"""
function Convergence(init;smooth=false,name=nothing,T=0.5)

	problems=Problem[] 	# array of problems 
	E0=Float64[]		# array or relative L^2 errors
	E1=Float64[]		# araay of relative H^1 errors
	
<<<<<<< Updated upstream
	if init == 1
		reference_problem=IntegrateSV(init=init,α=1.5,ϵ=1,L=π,N=2^15,T=0.05,dt = 1e-5,dealias=true,smooth=false,Ns=100,label="reference")
		Uref=reference_problem.data.U[end]/2^15

		for n=6:14
			p = IntegrateSV(init=init,α=1.5,ϵ=1,L=π,N=2^n,T=0.05,dt = 1e-5,dealias=true,smooth=smooth,Ns=1,label="N=2^$n")
=======
	if init == 1 # initial data is a heap of water
		# Solve the reference problem
		reference_problem=IntegrateSV(init=init,α=1.5,ϵ=1,L=π,M=2^15,T=T,dt = 1e-5,dealias=true,smooth=false,Ns=100,label="reference")
		# Save the (Fourier modes of) last computed solution
		Uref=reference_problem.data.U[end]/2^15

		for n=14:-1:6
			# Solve the problem with fewer Fourier modes
			p = IntegrateSV(init=init,α=1.5,ϵ=1,L=π,M=2^n,T=T,dt = 1e-5,dealias=true,smooth=smooth,Ns=1,label="2M=2^$n")
>>>>>>> Stashed changes
			push!(problems,p)
			# Save the (Fourier modes of) last computed solution
			U=[p.data.U[end][1] ; p.data.U[end][2]]/2^n
			# Extract the corresponding modes of the reference solution
			Nmodes=2^(n-1)
			Ucomp=[ Uref[1][1:Nmodes];  Uref[1][end-Nmodes+1:end] ;
						Uref[2][1:Nmodes];  Uref[2][end-Nmodes+1:end]  ]

			# Compute the relative errors
			push!(E0,norm(U-Ucomp,2)/norm(Ucomp,2))
			k=[0:Nmodes-1; -Nmodes:-1]
			K=[sqrt.(1 .+ k.^2);sqrt.(1 .+ k.^2)]
			push!(E1,norm(K.*U-K.*Ucomp,2)/norm(K.*Ucomp,2))
		end
	elseif init == 2

		for n=14:-1:6
			# Solve the reference problem
			reference_problem=IntegrateSV(init=init,N=(2^n÷3),h₀=0.5,v₀=2,ϵ=1,L=π,M=2^15,T=T,dt = 1e-5,dealias=true,smooth=smooth,Ns=1,label="reference")
			# Solve the problem with fewer Fourier modes
			p = IntegrateSV(init=init,h₀=0.5,v₀=2,ϵ=1,L=π,M=2^n,T=T,dt = 1e-5,dealias=true,smooth=smooth,Ns=1,label="2M=2^$n")

			# Save the (Fourier modes of) last computed solutions
			Uref=reference_problem.data.U[end]/2^15
			U=[p.data.U[end][1] ; p.data.U[end][2]]/2^n
			# Extract the corresponding modes of the reference solution
			Nmodes=2^(n-1)
			Ucomp=[ Uref[1][1:Nmodes];  Uref[1][end-Nmodes+1:end] ;
						Uref[2][1:Nmodes];  Uref[2][end-Nmodes+1:end]  ]

			# Compute the relative errors
			push!(E0,norm(U-Ucomp,2)/norm(Ucomp,2))
			k=[0:Nmodes-1; -Nmodes:-1]
			K=[sqrt.(1 .+ k.^2);sqrt.(1 .+ k.^2)]
			push!(E1,norm(K.*U-K.*Ucomp,2)/norm(K.*Ucomp,2))
		end
	end
	return reference_problem,E0,E1

end


#--- Figures

#------- Experiments with a heap of water
function Figure1()
# Plot the initial data (heap of water) and the decay of Fourier coefficients
# Build the problems, the errors measure the decay of Fourier coefficients
reference_problem,E0,E1=Convergence(1;smooth=false,T=0)	

# Extract initial data of the reference problem
η0,v0,x=reference_problem.model.mapfro(reference_problem.data.U[1])
# Plot initial data
Fig1a=plot(x,[η0 v0],label=["\$\\eta^0\$" "\$u^0\$"],xlabel="\$x\$")
savefig(Fig1a,"Fig1a.pdf");savefig(Fig1a,"Fig1a.svg");
# Plot convergence rates
<<<<<<< Updated upstream
Fig1b=plot(;xlabel="\$N\$",axis=:log)
Ns=2 .^(6:14);
scatter!(Fig1b,Ns,E1,label="\$E_1\$",color=1)
scatter!(Fig1b,Ns,E0,label="\$E_0\$",color=2)
plot!(Fig1b,Ns,Ns.^(-1),label="",color=1)
plot!(Fig1b,Ns,Ns.^(-2),label="",color=2)
=======
Fig1b=plot(;xlabel="\$2M\$",axis=:log)
n=14:-1:6;M=2 .^n;
scatter!(Fig1b,M,E1,label="\$E_1\$",color=1)
scatter!(Fig1b,M,E0,label="\$E_0\$",color=2)
plot!(Fig1b,M,M.^(-1),label="",color=1)
plot!(Fig1b,M,M.^(-2),label="",color=2)
>>>>>>> Stashed changes
savefig(Fig1b,"Fig1b.pdf");savefig(Fig1b,"Fig1b.svg");
end

#--- Experiment 1 : heap of water, sharp low-pass filter
function Figure2()
# Solve the problems, save the errors and the reference solution
reference_problem,E0,E1=Convergence(1;smooth=false)	

#--- Experiment 2 : heap of water, smooth low-pass filter
# Solve the problems, save the errors and the reference solution
reference_problem,E0_smooth,E1_smooth=Convergence(1;smooth=true)	

# Plot convergence rates
<<<<<<< Updated upstream
Fig2=plot(;xlabel="\$N\$",axis=:log)
Ns=2 .^(6:14)
scatter!(Fig2,Ns,E1,label="\$E_1\$, sharp low-pass filter",color=1)
scatter!(Fig2,Ns,E0,label="\$E_0\$, sharp low-pass filter",color=2)
scatter!(Fig2,Ns,E1_smooth,label="\$E_1\$, smooth low-pass filter",color=3)
scatter!(Fig2,Ns,E0_smooth,label="\$E_0\$, smooth low-pass filter",color=4)
=======
Fig2=plot(;xlabel="\$2M\$",axis=:log)
n=14:-1:6;M=2 .^n;
scatter!(Fig2,M,E1,label="\$E_1\$, sharp low-pass filter",color=1)
scatter!(Fig2,M,E0,label="\$E_0\$, sharp low-pass filter",color=2)
scatter!(Fig2,M,E1_smooth,label="\$E_1\$, smooth low-pass filter",color=3)
scatter!(Fig2,M,E0_smooth,label="\$E_0\$, smooth low-pass filter",color=4)
>>>>>>> Stashed changes

plot!(Fig2,M,M.^(-1),label="",color=1)
plot!(Fig2,M,M.^(-2),label="",color=2)

savefig(Fig2,"Fig2.pdf");savefig(Fig1b,"Fig2.svg");

# Table : experimental order of Convergence
return round.(log2.(E0[1:end-1]./E0[2:end]),digits=2),
		round.(log2.(E1[1:end-1]./E1[2:end]),digits=2),
		round.(log2.(E0_smooth[1:end-1]./E0_smooth[2:end]),digits=2),
		round.(log2.(E1_smooth[1:end-1]./E1_smooth[2:end]),digits=2)
end


#--------- Experiments with a high frequency component
function Figure3()
# Plot the initial data (heap of water) and the decay of Fourier coefficients
# Build the problems, the errors measure the decay of Fourier coefficients
reference_problem,E0,E1=Convergence(2;smooth=false,T=0)	

# Extract initial data of the first problem
η0,v0,x=reference_problem.model.mapfro(reference_problem.data.U[1])

# Plot initial data
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

#--- Experiment 3 : high frequency component, sharp low-pass filter
# Solve the problems, save the errors and the reference solution
reference_problem,E0,E1=Convergence(2;smooth=false)	

<<<<<<< Updated upstream
# Plot convergence rates
Fig3b=plot(;xlabel="\$N\$",axis=:log)
=======
#--- Experiment 4 : high frequency component, smooth low-pass filter
reference_problem,E0_smooth,E1_smooth=Convergence(2;smooth=true)	

# Plot convergence rates
Fig3b=plot(;xlabel="\$2M\$",axis=:log,legend=:bottomleft)
>>>>>>> Stashed changes
Ns=2 .^(6:14);
scatter!(Fig3b,Ns,E1[end:-1:1],label="\$E_1\$",color=1)
scatter!(Fig3b,Ns,E0[end:-1:1],label="\$E_0\$",color=2)
plot!(Fig3b,Ns,Ns.^(-1),label="",color=1)
plot!(Fig3b,Ns,Ns.^(-2),label="",color=2)
savefig(Fig3b,"Fig3b.pdf");savefig(Fig1b,"Fig3b.svg");
end

<<<<<<< Updated upstream

# Some experiments to play around experiment 3: sharp cut-off vs a reference solution
n=10
p = IntegrateSV(init=2,h₀=0.5,v₀=2,ϵ=1,L=π,N=2^n,T=0.05,dt = 1e-5,dealias=true,smooth=false,Ns=1,label="N=2^$n")
reference_problem=IntegrateSV(init=2,M=(2^n÷3),h₀=0.5,v₀=2,ϵ=1,L=π,N=2^12,T=0.05,dt = 1e-5,dealias=true,smooth=false,Ns=1,label="N=2^15")
plot(p,T=0,var=:fourier_velocity)
plot!(reference_problem,T=0,var=:fourier_velocity)
xlims!(-5,5)
η,v,x=solution(p;T=0)
ηref,vref,xref=solution(reference_problem;T=0)
plot(x,v.-2*sin.(x))
plot!(xref,vref.-2*sin.(xref))

plot(p;var=[:surface :velocity :fourier :fourier_velocity])
plot!(reference_problem;var=[:surface :velocity :fourier :fourier_velocity])
=======
#--- Experiment 5 : instability in zero-depth situation
function Figure4()
# Solve problem with vanishing depth and smooth cut-off and M=2^10 modes
pb_smooth_M10=IntegrateSV(init=2,h₀=0,v₀=2,ϵ=1,L=π,M=2^10,T=0.1,dt = 1e-5,dealias=true,smooth=true,Ns=1,label="smooth")
# Solve problem with vanishing depth and sharp cut-off and M=2^10 modes
pb_sharp_M10=IntegrateSV(init=2,h₀=0,v₀=2,ϵ=1,L=π,M=2^10,T=0.1,dt = 1e-5,dealias=true,smooth=false,Ns=1,label="sharp")
>>>>>>> Stashed changes

# Plot solutions
plot(pb_sharp_M10,var=[:surface :velocity :fourier_velocity])
plot!(pb_smooth_M10,var=[:surface :velocity :fourier_velocity])

# Plot second derivatives of solutions
η_M10,v_M10,x=solution(pb_sharp_M10)
ηs_M10,vs_M10,=solution(pb_smooth_M10)
∂=1im*Mesh(x).k
diff(η;n=1)=real.(ifft(∂.^n.*fft(η)))

Fig4a=plot(x,diff(v_M10;n=2),label="sharp")
plot!(x,diff(vs_M10;n=2),label="smooth")
title!("\$\\partial_x^2u_N, \\quad    2M = 2^{10}\$")
xlabel!("\$x\$")
savefig(Fig4a,"Fig4a.pdf");savefig(Fig4a,"Fig4a.svg");



# Solve problem with vanishing depth and smooth cut-off and M=2^12 modes
pb_smooth_M12=IntegrateSV(init=2,h₀=0,v₀=2,ϵ=1,L=π,M=2^12,T=0.1,dt = 1e-5,dealias=true,smooth=true,Ns=1,label="smooth")
# Solve problem with vanishing depth and sharp cut-off and M=2^12 modes
pb_sharp_M12=IntegrateSV(init=2,h₀=0,v₀=2,ϵ=1,L=π,M=2^12,T=0.1,dt = 1e-5,dealias=true,smooth=false,Ns=1,label="sharp")

# Plot solutions
plot(pb_sharp_M12,var=[:surface :velocity :fourier_velocity])
plot!(pb_smooth_M12,var=[:surface :velocity :fourier_velocity])

# Plot second derivatives of solutions
η_M12,v_M12,x=solution(pb_sharp_M12)
ηs_M12,vs_M12,=solution(pb_smooth_M12)
∂=1im*Mesh(x).k
diff(η;n=1)=real.(ifft(∂.^n.*fft(η)))

<<<<<<< Updated upstream
Fig4b=plot(x,diff(v;n=2),label="sharp")
ηs,vs,=solution(pb_smooth)
plot!(x,diff(vs;n=2),label="smooth")
title!("\$\\partial_x^2v\$")
xlabel!("\$x\$")
savefig(Fig4b,"Fig4b.pdf");savefig(Fig4b,"Fig4b.svg");
=======
Fig4b=plot(x,diff(v_M12;n=2),label="sharp")
plot!(x,diff(vs_M12;n=2),label="smooth")
title!("\$\\partial_x^2u_N, \\quad    2M = 2^{12}\$")
xlabel!("\$x\$")
savefig(Fig4b,"Fig4b.pdf");savefig(Fig4b,"Fig4b.svg");
end

#--- Experiment 6 (not shown in the manuscript) : simluations after wavebreaking
function FigureNone()
# Solve problem with smooth simple wave initial data and M=2^10 modes
pb_sharp_M10=IntegrateSV(init=3,ϵ=1,L=π,M=2^10,T=1,dt = 1e-5,dealias=true,smooth=false,Ns=100,label="sharp")
pb_smooth_M10=IntegrateSV(init=3,ϵ=1,L=π,M=2^10,T=1,dt = 1e-5,dealias=true,smooth=true,Ns=100,label="smooth")
plot([pb_sharp_M10,pb_smooth_M10],T=1)

# Solve problem with smooth simple wave initial data and M=2^12 modes
pb_sharp_M12=IntegrateSV(init=3,ϵ=1,L=π,M=2^12,T=1,dt = 1e-5,dealias=true,smooth=false,Ns=100,label="sharp")
pb_smooth_M12=IntegrateSV(init=3,ϵ=1,L=π,M=2^12,T=1,dt = 1e-5,dealias=true,smooth=true,Ns=100,label="smooth")
plot([pb_sharp_M12,pb_smooth_M12],T=1)
end


#------------------- 
#--- 2D PROBLEMS ---
#-------------------

"""
	IntegrateSV2D(;init,args)

Integrate in time two-dimensional Saint-Venant (shallow water) system.

All arguments are optional:
- `(N,s,u₀,v₀,u₁,v₁) describe parameters of the initial initial data (by default, `N` is 2/3 of the number of collocation points),
- `ϵ` the nonlinearity parameter (default is `1`),
- `L` the half-length of the mesh (default is `π`),
- `M` the number of collocation points (default is `2^9`),
- `T` the final time of integration (default is `0.1`),
- `dt` the timestep (default is `5e-5`),
- `dealias`: no dealisasing if set to `0` or `false`, otherwise `1/(3*dealias)` modes are set to `0` (corresponding to standard 2/3 Orszag rule if `dealias` is set to `1` or `true`, which is default),
- `smooth`: A smooth low-pass filter (whose scaling is defined by `N`) if set to `0` or `false` (default), otherwise only `2/(3*dealias)*(1-smooth/2)` modes are kept untouched,
- `hamiltonian`: if `true` (default is `false`) then the Hamiltonian Saint-Venant system is computed,
- `Ns` the number of stored computed times (default is only initial and final times),
- `label`: a label (string) for future use (default is ""Saint-Venant"").


Return `problem` of the solved problem 
"""
function IntegrateSV2D(;N=nothing,s=2,h₀=1/2,u₀=1/2,v₀=-1/2,u₁=1,v₁=1,ϵ=1,L=π,M=2^9,T=0.1,dt =5e-5,dealias=true,smooth=false,hamiltonian=false,Ns=1,label="Saint-Venant")
	if isnothing(Ns)
		param = ( ϵ  = ϵ,
				N  = M, L  = L,
	            T  = T, dt = dt )
	else
		param = ( ϵ  = ϵ,
				N  = M, L  = L,
	            T  = T, dt = dt,
				Ns = Ns )
	end
>>>>>>> Stashed changes

	mesh=Mesh(param)
	if isnothing(N) N=(M÷3) end
	
	ζ(x,y) = (h₀-1) * cos.(x).*cos.(y)'; 
    ux(x,y) = u₀* cos.(y').*sin.(x) .+ u₁*cos.(N*y').*sin.(N*x)/N^s; 
    uy(x,y) = v₀* sin.(y').*cos.(x) .+ v₁*sin.(N*y').*cos.(N*x)/N^s;

	init = Init2D(ζ, ux, uy);

	model=SaintVenant2D_fast(param; dealias = dealias , hamiltonian = hamiltonian, smooth=smooth)
    solver=RK4(model.mapto(init))
    problem = Problem(  model, init, param ; solver=solver )
	
	solve!( problem )

<<<<<<< Updated upstream
#--- Experiment 5 : instability ?
pb_ref=IntegrateSV(init=2,h₀=0.5,v₀=2,ϵ=1,L=π,N=2^13,T=0.1,dt = 1e-5,dealias=true,smooth=false,Ns=1,label="ref")
pb_ins=IntegrateSV(init=3,h₀=0.5,v₀=2,ϵ=1,L=π,N=2^13,T=0.1,dt = 1e-5,dealias=true,smooth=false,Ns=1,label="ins")
Fig5a=plot(pb_ref,var=[:surface :velocity :fourier_velocity])
plot!(pb_ins,var=[:surface :velocity :fourier_velocity])
#savefig(Fig4a,"Fig4a.pdf");savefig(Fig4a,"Fig4a.svg");
=======
	# Check energy preservation
	e(η,vx,vy) = 1/2* ( η.^2 .+ (1 .+ϵ*η).*(vx.^2 .+vy.^2) );
	e(η1,vx1,vy1,η2,vx2,vy2) = 1/2* ( (η1-η2).*(η1+η2) .+ (1 .+ϵ*η2).*(vx1-vx2).*(vx1+vx2) .+ (1 .+ϵ*η2).*(vy1-vy2).*(vy1+vy2)
						.+ϵ*(η1-η2).*vx1.^2 .+ ϵ*(η1-η2).*vy1.^2 );
	
	(ηfin,vxfin,vyfin)   =  solution(problem)
	(η0,vx0,vy0)=(init.η(mesh.x,mesh.x),init.vx(mesh.x,mesh.x),init.vy(mesh.x,mesh.x))

	error_energy = abs(sum(e(ηfin,vxfin,vyfin,η0,vx0,vy0))/sum(e(η0,vx0,vy0)))
	@info "normalized preservation of the total energy: $error_energy\n"

	return problem
end


"""
	Convergence(init;smooth,T=0.5)

Convergence rate associated with numerical experiments.
The Saint-Venant (shallow water) system is numerically integrated with several values for the number of collocation points.

If optional argument `smooth` is `true` (default is `false`): a smooth low-pass filter is used.
If optional argument `hamiltonian` is `true` (default is `false`): the Hamiltonian Saint-Venant system is computed,

The optional argument `T` is the final time of computation (default is `0.1`).

Return the reference problem and relative errors.
"""
function Convergence2D(;smooth=false,hamiltonian=false,name=nothing,T=0.1)

	problems=Problem[] 	# array of problems 
	E0=Float64[]		# array or relative L^2 errors
	E1=Float64[]		# araay of relative H^1 errors
	
	

	for n=9:-1:6
		# Solve the reference problem
		reference_problem=IntegrateSV2D(M=2^10,N = (2^n)÷3,T=T,smooth=smooth,hamiltonian=hamiltonian,Ns=1,label="reference")
		# Save the (Fourier modes of) last computed solution
		Uref=reference_problem.data.U[end]/2^20
        Uref1 = fftshift(Uref[1])
        Uref2 = fftshift(Uref[2])
        Uref3 = fftshift(Uref[3])

		# Solve the problem with fewer Fourier modes
		p = IntegrateSV2D(M=2^n,N = (2^n)÷3, T=T,smooth=smooth,hamiltonian=hamiltonian,Ns=1,label="2M=2^$n")
        push!(problems,p)
   		# Save the (Fourier modes of) last computed solution
		U = p.data.U[end]/2^(2n)
        U1 = fftshift(U[1])
        U2 = fftshift(U[2])
        U3 = fftshift(U[3])
        # Extract the corresponding modes of the reference solution
		Nmodes=2^(n-1);Nref=2^9;
		Ucomp1 = Uref1[Nref-Nmodes+1:Nref+Nmodes, Nref-Nmodes+1:Nref+Nmodes];
        Ucomp2 = Uref2[Nref-Nmodes+1:Nref+Nmodes, Nref-Nmodes+1:Nref+Nmodes];
        Ucomp3 = Uref3[Nref-Nmodes+1:Nref+Nmodes, Nref-Nmodes+1:Nref+Nmodes];

		# Compute the relative errors
		kx = range(-Nmodes,Nmodes-1)
        ky = kx'
        k = sqrt.(1 .+ kx.^2 .+ ky.^2)
        
        push!(E0,norm([norm(U1-Ucomp1,2);norm(U2-Ucomp2,2);norm(U3 - Ucomp3,2)],2)/norm([norm(Ucomp1,2);norm(Ucomp2,2);norm(Ucomp3,2)],2))
        push!(E1,norm([norm(k.*U1-k.*Ucomp1,2);norm(k.*U2-k.*Ucomp2,2);norm(k.*U3 - k.*Ucomp3,2)],2)/norm([norm(k.*Ucomp1,2);norm(k.*Ucomp2,2);norm(k.*Ucomp3,2)],2))
	end

	return E0,E1

end

#--- Figures

#------- Experiment 1. Convergence rate for initial data satisfying the strong hyperbolicity condition
function Figure5()
#--- Non-Hamiltonian system

# Solve the problems, save the errors
E0,E1=Convergence2D(;T=0.1,smooth=false, hamiltonian=false)
E0s,E1s=Convergence2D(;T=0.1,smooth=true, hamiltonian=false)

# Plot convergence rates
Fig5a=plot(;xlabel="\$M\$",axis=:log,legend=:bottomleft)
n=9:-1:6;Ns=2 .^n;
scatter!(Fig5a,Ns,E1,label="\$E_1\$, sharp low-pass filter",color=1)
scatter!(Fig5a,Ns,E0,label="\$E_0\$, sharp low-pass filter",color=2)
scatter!(Fig5a,Ns,E1s,label="\$E_1\$, smooth low-pass filter",color=3)
scatter!(Fig5a,Ns,E0s,label="\$E_0\$, smooth low-pass filter",color=4)
plot!(Fig5a,Ns,4 .* Ns.^(-1),label="",color=1)
plot!(Fig5a,Ns,12 .*Ns.^(-2),label="",color=2)

savefig(Fig5a,"Fig5a.pdf");savefig(Fig5a,"Fig5a.svg");

>>>>>>> Stashed changes


#--- Hamiltonian system

# Solve the problems, save the errors
E0,E1=Convergence2D(;T=0.1,smooth=false, hamiltonian=true)
E0s,E1s=Convergence2D(;T=0.1,smooth=true, hamiltonian=true)

# Plot convergence rates
Fig5b=plot(;xlabel="\$M\$",axis=:log,legend=:bottomleft)
n=9:-1:6;Ns=2 .^n;
scatter!(Fig5b,Ns,E1,label="\$E_1\$, sharp low-pass filter",color=1)
scatter!(Fig5b,Ns,E0,label="\$E_0\$, sharp low-pass filter",color=2)
scatter!(Fig5b,Ns,E1s,label="\$E_1\$, smooth low-pass filter",color=3)
scatter!(Fig5b,Ns,E0s,label="\$E_0\$, smooth low-pass filter",color=4)
plot!(Fig5b,Ns,4 .* Ns.^(-1),label="",color=1)
plot!(Fig5b,Ns,12 .*Ns.^(-2),label="",color=2)

savefig(Fig5b,"Fig5b.pdf");savefig(Fig5b,"Fig5b.svg");
end

#------- Experiment 2. Negative-depth 
function Figure6()
##---- Solve problems (M=2^9)
problem = IntegrateSV2D(;h₀=-0.1,u₀=1/2,v₀=-1/2,u₁=1,v₁=1,M=2^9,T=0.1,dt =5e-5,smooth=false,hamiltonian=false,Ns=1)
problem_smooth =  IntegrateSV2D(;h₀=-0.1,u₀=1/2,v₀=-1/2,u₁=1,v₁=1,M=2^9,T=0.1,dt =5e-5,smooth=true,hamiltonian=false,Ns=1)

η,vx,vy,x,y = solution(problem)
ηs,vxs,vys, = solution(problem_smooth)

mesh=Mesh(x)
kx=mesh.k;ky=mesh.k;x=mesh.x;y=mesh.x;
∂y = 1im * ky' ;
∂x = 1im * kx  ;

function ∂fx( f )
    real.(ifft(∂x .* fft(f, 1), 1))
end
function ∂fy( f )
    real.(ifft(∂y .* fft(f, 2), 2))
end

# sharp low-pass filter, non-hamiltonian system
Fig6a=plot(x,y,∂fx(∂fx(vx)),st=:surface,label="sharp")
plot!(xlabel="x", ylabel="y")
title!("\$\\partial_x^2u_N, \\quad    2M = 2^{9}\$")
savefig(Fig6a,"Fig6a.pdf");savefig(Fig6a,"Fig6a.svg");

# smooth low-pass filter, non-hamiltonian system
Fig6b=plot(x,y,∂fx(∂fx(vxs)),st=:surface,label="smooth")
plot!(xlabel="x", ylabel="y")
title!("\$\\partial_x^2u_N, \\quad    2M = 2^{9}\$")
savefig(Fig6b,"Fig6b.pdf");savefig(Fig6b,"Fig6b.svg");

##---- Solve problems (M=2^10)
problem = IntegrateSV2D(;h₀=-0.1,u₀=1/2,v₀=-1/2,u₁=1,v₁=1,M=2^10,T=0.1,dt =5e-5,smooth=false,hamiltonian=false,Ns=1)
problem_smooth =  IntegrateSV2D(;h₀=-0.1,u₀=1/2,v₀=-1/2,u₁=1,v₁=1,M=2^10,T=0.1,dt =5e-5,smooth=true,hamiltonian=false,Ns=1)

η,vx,vy,x,y = solution(problem)
ηs,vxs,vys, = solution(problem_smooth)

mesh=Mesh(x)
kx=mesh.k;ky=mesh.k;x=mesh.x;y=mesh.x;
∂y = 1im * ky' ;
∂x = 1im * kx  ;

function ∂fx( f )
    real.(ifft(∂x .* fft(f, 1), 1))
end
function ∂fy( f )
    real.(ifft(∂y .* fft(f, 2), 2))
end

# sharp low-pass filter, non-hamiltonian system
Fig6c=plot(x,y,∂fx(∂fx(vx)),st=:surface,label="sharp")
plot!(xlabel="x", ylabel="y")
title!("\$\\partial_x^2u_N, \\quad    2M = 2^{10}\$")
savefig(Fig6c,"Fig6c.pdf");savefig(Fig6c,"Fig6c.svg");

# smooth low-pass filter, non-hamiltonian system
Fig6d=plot(x,y,∂fx(∂fx(vxs)),st=:surface,label="smooth")
plot!(xlabel="x", ylabel="y")
title!("\$\\partial_x^2u_N, \\quad    2M = 2^{10}\$")
savefig(Fig6d,"Fig6d.pdf");savefig(Fig6d,"Fig6d.svg");
end

#------- Experiment 3. Hamiltonian system 
function Figure7()
##---- Solve problems (M=2^9)
problem = IntegrateSV2D(;h₀=0.5,u₀=2,v₀=-2,u₁=1,v₁=-1,M=2^9,T=0.1,dt =5e-5,smooth=false,hamiltonian=true,Ns=1)
problem_smooth =  IntegrateSV2D(;h₀=0.5,u₀=2,v₀=-2,u₁=1,v₁=-1,M=2^9,T=0.1,dt =5e-5,smooth=true,hamiltonian=true,Ns=1)

η,vx,vy,x,y = solution(problem)
ηs,vxs,vys, = solution(problem_smooth)

mesh=Mesh(x)
kx=mesh.k;ky=mesh.k;x=mesh.x;y=mesh.x;
∂y = 1im * ky' ;
∂x = 1im * kx  ;

function ∂fx( f )
    real.(ifft(∂x .* fft(f, 1), 1))
end
function ∂fy( f )
    real.(ifft(∂y .* fft(f, 2), 2))
end

# sharp low-pass filter, Hamiltonian system
Fig7a=plot(x,y,∂fx(∂fx(vx)),st=:surface,label="sharp")
plot!(xlabel="x", ylabel="y")
title!("\$\\partial_x^2u_N, \\quad    2M = 2^{9}\$")
savefig(Fig7a,"Fig7a.pdf");savefig(Fig7a,"Fig7a.svg");

# smooth low-pass filter, non-hamiltonian system
Fig7b=plot(x,y,∂fx(∂fx(vxs)),st=:surface,label="smooth")
plot!(xlabel="x", ylabel="y")
title!("\$\\partial_x^2u_N, \\quad    2M = 2^{9}\$")
savefig(Fig7b,"Fig7b.pdf");savefig(Fig7b,"Fig7b.svg");

##---- Solve problems (M=2^10)
problem = IntegrateSV2D(;h₀=0.5,u₀=2,v₀=-2,u₁=1,v₁=-1,M=2^10,T=0.1,dt =5e-5,smooth=false,hamiltonian=true,Ns=1)
problem_smooth =  IntegrateSV2D(;h₀=0.5,u₀=2,v₀=-2,u₁=1,v₁=-1,M=2^10,T=0.1,dt =5e-5,smooth=true,hamiltonian=true,Ns=1)

η,vx,vy,x,y = solution(problem)
ηs,vxs,vys, = solution(problem_smooth)

mesh=Mesh(x)
kx=mesh.k;ky=mesh.k;x=mesh.x;y=mesh.x;
∂y = 1im * ky' ;
∂x = 1im * kx  ;

function ∂fx( f )
    real.(ifft(∂x .* fft(f, 1), 1))
end
function ∂fy( f )
    real.(ifft(∂y .* fft(f, 2), 2))
end

# sharp low-pass filter, non-hamiltonian system
Fig7c=plot(x,y,∂fx(∂fx(vx)),st=:surface,label="sharp")
plot!(xlabel="x", ylabel="y")
title!("\$\\partial_x^2u_N, \\quad    2M = 2^{10}\$")
savefig(Fig7c,"Fig7c.pdf");savefig(Fig7c,"Fig7c.svg");

# smooth low-pass filter, non-hamiltonian system
Fig7d=plot(x,y,∂fx(∂fx(vxs)),st=:surface,label="smooth")
plot!(xlabel="x", ylabel="y")
title!("\$\\partial_x^2u_N, \\quad    2M = 2^{10}\$")
savefig(Fig7d,"Fig7d.pdf");savefig(Fig7d,"Fig7d.svg");
end
nothing