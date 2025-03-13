# #
# Reproduces the figures in
# [V. Duchêne and B. Melinand](https://doi.org/10.2140/paa.2024.6.73)
# on the quadratic pseudo-spectral method (WW2)
# #
export IntegrateWW2,Figure
using WaterWaves1D,FFTW,Plots,LinearAlgebra,ProgressMeter;
#using JLD #(uncomment if using @save)

#--- Integration

"""
	IntegrateWW2(;init,args)

Integrate in time the WW2 system with an initial data depending on the provided `init`
- if `init=1`, then surface deformation `η(t=0,x)=exp(-x^p)` and velocity `v(t=0,x)=0` (with `p` provided as an optional argument, by default `p=2`)
- if `init=2`, then surface deformation `η(t=0,x)=exp(-x^2)` and velocity `v(t=0,x)=exp.(-x^2)*(sin(x)+sin(K*x)/K^2)` (with `K` provided as an optional argument, by default `K=100`)

Other arguments are optional:
- `μ` the shallowness parameter (default is `1`),
- `ϵ` the nonlinearity parameter (default is `0.1`),
- `L` the half-length of the mesh (default is `20`),
- `N` the number of collocation points (default is `2^10`),
- `T` the final time of integration (default is `10`),
- `dt` the timestep (default is `0.001`),
- `dealias`: dealiasing with Orszag rule `1-dealias/(dealias+2)` (default is `1`, i.e. 2/3 rule, `0` means no dealiasing);
- `δ` the strength of the rectifier (default is `0.001`),
- `m` the order of the rectifier (as a regularizing operator, default is `-1`),
- `Ns` the number of stored computed times (default is all times).


Return `(problem,blowup_time,blowup,error_energy)` where
- `problem` contains all the information,
- `blowup_time` is the first time with NaN values (final time if there is none)
- `blowup` is a boolean indicating if NaN values occured
- `error_energy` is the relative energy preservation between first and final time
"""
function IntegrateWW2(;init=1,μ=1,ϵ=0.1,L=20,N=2^10,T=10,dt = 0.001,dealias=1,δ=0.001,m=-1,K=100,p=2,Ns=nothing)
	if isnothing(Ns)
		param = ( μ  = μ, ϵ  = ϵ,
				N  = N, L  = L,
	            T  = T, dt = dt )
	else
		param = ( μ  = μ, ϵ  = ϵ,
				N  = N, L  = L,
	            T  = T, dt = dt,
				Ns = Ns )
	end

	mesh=Mesh(param)
	if init == 1
		init = Init(x->exp.(-abs.(x).^p),x->zero(x))
	elseif init == 2
		init = Init(x->zero(x), x->exp.(-x.^2).*(sin.(x).+sin.(K*x)/K^2))
	else
		@error "argument init must be 1 or 2"
	end

	model = WWn(param;n = 2, δ = δ, m = m, dealias = dealias)
	problem = Problem(model, init, param)
	solve!( problem )

	(ηfin,vfin)   =  solution(problem)

	# check energy preservation
	x=mesh.x;k=mesh.k;
	Tmu = -1im*tanh.(sqrt(μ)*k)
	tmu = tanh.(sqrt(μ)*k)./k;tmu[1]=sqrt(μ);
	∂ₓ= 1im*k


	e(η,v) = 1/2* ( η.^2 .+ v.*ifft(tmu.*fft(v)) +ϵ*η.*(v.^2-ifft(Tmu.*fft(v)).^2) );
	e(η1,v1,η2,v2) = 1/2* ( (η1-η2).*(η1+η2) .+ (v1-v2).*ifft(tmu.*fft(v1+v2))
						+ϵ*(η1-η2).*(v1.^2-ifft(Tmu.*fft(v1)).^2)
						+ϵ*η2.*((v1-v2).*(v1+v2)-ifft(Tmu.*fft(v1-v2)).*ifft(Tmu.*fft(v1+v2))) );
	η0=init.η(x);v0=init.v(x)
	#error_energy_1 = abs(sum(e(ηfin,vfin)-e(η0,v0))/sum(e(η0,v0)))
	error_energy = abs(sum(e(ηfin,vfin,η0,v0))/sum(e(η0,v0)))
	@info "normalized preservation of the total energy: $error_energy\n"

	# compute blowup time (if any)
	blowup_index=0;blowup=false;
	for j in 1:problem.times.Ns
		blowup_index+=1
		if isnan(problem.data.U[j][1][1])
			blowup=true
			break
		end
	end
	blowup_time=	problem.times.ts[blowup_index]
	if blowup_time<T
		@info "The solution blowed at time t=$blowup_time (first NaN occurence)\n"
	end
	#display(plt)
	return problem,blowup_time,blowup,error_energy
end

#--- Figures

"""
	Figure(scenario;name,anim)

Several numerical experiments,
depending on the argument `scenario` (between 1 and 12)
corresponding to different figures in [V. Duchêne and B. Melinand](https://doi.org/10.2140/paa.2024.6.73)

- `scenario∈[1,7]`: spurious instability formation for smooth initial data.
    - `scenario=1`: no instability even without dealiasing and rectification, if small number of modes.
    - `scenario=2`: instabilities without dealiasing and rectification, if larger number of modes.
    - `scenario=3`: no instability with dealiasing and without rectification, if small number of modes.
    - `scenario=4`: instabilities with dealiasing and without rectification, if larger number of modes.
    - `scenario=5`: no instability in the presence of sufficiently regularizing rectifiers.
    - `scenario=6`: instabilities in the presence of insufficiently regularizing rectifiers.
    - `scenario=7`: mild instabilities in the presence of insufficiently strong rectifiers (δ small).
- `scenario∈[8,10]`: blowups for an initial data with a low-frequency and a high-frequency (with wavenumber K) component.
    - `scenario=8`: dependency of blowup times with respect to K (also `scenario=8.5` with a cut-off proportional to K)
    - `scenario=9`: dependency of blowup times with respect to ϵ.
    - `scenario=10`: an example of blowup with two different dealiasing.
- `scenario=11`: critical rectifier strength (δ) as a function of ϵ.
- `scenario≈12`: error of the rectified model (WW2) with respect to the water waves system, depending on δ and ϵ.
    - `scenario=12.1`: with initial data exp(-|x|)
    - `scenario=12.2`: with initial data exp(-|x|^2)
    - `scenario=12.3`: with initial data exp(-|x|^3)


Optional arguments are
- `compression`, an integer `m` used to reduce the size of figures by plotting only one every `m` points (when relevant),
- `name` (a string) used to save figures,
- `anim` generate animations (when relevant) if set to `true`.


Return relevant plots and problems depending on the situation.
"""
function Figure(scenario;compression=false,name=nothing,anim=false)

	if scenario == 1
		problem,=IntegrateWW2(init=1,μ=1,ϵ=0.1,L=20,N=2^9,T=10,dt = 0.001,dealias=0,δ=0,m=-1)
		plt=plot(problem;var=[:surface,:fourier],compression=compression,label="")
		if !isnothing(name)
			savefig(plt,"$name.pdf");savefig(plt,"$name.svg");
			if anim
				anim = @animate for t in LinRange(0,problem.times.tfin,101)
					plt = plot(problem;var=[:surface,:fourier],T=t,compression=compression)
					ylims!(plt[1],(-0.5,1.1))
				end
				gif(anim, "$name.gif", fps = 15)
			end
 		end
		return problem,plt

	elseif scenario == 2
		problem,=IntegrateWW2(init=1,μ=1,ϵ=0.1,L=20,N=2^11,T=1.5,dt = 0.001,dealias=0,δ=0,m=-1)
		plt=plot(problem,var=[:surface,:fourier],T=1.2;compression=compression,label="")
		if !isnothing(name)
			savefig(plt,"$name.pdf");savefig(plt,"$name.svg");
			if anim
				anim = @animate for t in LinRange(0,problem.times.tfin,101)
					plt = plot(problem;var=[:surface,:fourier],T=t,compression=compression)
					ylims!(plt[1],(-0.5,1.1))
				end
				gif(anim, "$name.gif", fps = 15)			
			end
		end
		return problem,plt

	elseif scenario == 3
		problem,=IntegrateWW2(init=1,μ=1,ϵ=0.1,L=20,N=2^12,T=10,dt = 0.001,dealias=1,δ=0,m=-1)
		plt=plot(problem;var=[:surface,:fourier],compression=compression,label="")
		if !isnothing(name)
			savefig(plt,"$name.pdf");savefig(plt,"$name.svg");
			if anim
				anim = @animate for t in LinRange(0,problem.times.tfin,101)
					plt = plot(problem;var=[:surface,:fourier],T=t,compression=compression)
					ylims!(plt[1],(-0.5,1.1))
				end
				gif(anim, "$name.gif", fps = 15)			
			end
		end
		return problem,plt

	elseif scenario == 4
		problem,=IntegrateWW2(init=1,μ=1,ϵ=0.1,L=20,N=2^14,T=1.5,dt = 0.001,dealias=1,δ=0,m=-1)
		plt=plot(problem,var=[:surface,:fourier],T=1.3;compression=compression,label="")
		if !isnothing(name)
			savefig(plt,"$name.pdf");savefig(plt,"$name.svg");
			if anim
				anim = @animate for t in LinRange(0,problem.times.tfin,101)
					plt = plot(problem;var=[:surface,:fourier],T=t,compression=compression)
					ylims!(plt[1],(-0.5,1.1))
				end
				gif(anim, "$name.gif", fps = 15)			
			end
 		end
		return problem,plt

	elseif scenario == 5
		if anim Ns=nothing else Ns=1 end
		problem0,=IntegrateWW2(init=1,μ=1,ϵ=0.1,L=20,N=2^18,T=10,dt = 0.01,dealias=1,δ=0.01,m=-1/2,Ns=Ns)
		problem1,=IntegrateWW2(init=1,μ=1,ϵ=0.1,L=20,N=2^18,T=10,dt = 0.01,dealias=1,δ=0.01,m=-1,Ns=Ns)
		plt0=plot(problem0;var=[:surface,:fourier],compression=compression,label="")
		plt1=plot(problem1;var=[:surface,:fourier],compression=compression,label="")
		if !isnothing(name)
			savefig(plt0,"$name-r12.pdf");savefig(plt0,"$name-r12.svg");
			savefig(plt1,"$name-r1.pdf");savefig(plt1,"$name-r1.svg");
			if anim
				anim = @animate for t in LinRange(0,problem0.times.tfin,101)
					plt = plot([problem0,problem1];var=[:surface,:fourier],T=t,compression=compression,label=["m=-1/2" "m=-1"])
					ylims!(plt[1],(-0.5,1.1))
				end
				gif(anim, "$name.gif", fps = 15)
			end
		end
		return problem0,problem1,plt0,plt1

	elseif scenario == 6
		if anim Ns=nothing else Ns=10 end
		problem0,=IntegrateWW2(init=1,μ=1,ϵ=0.1,L=20,N=2^18,T=1,dt = 0.01,dealias=1,δ=0.01,m=-1/4,Ns=Ns)
		problem1,=IntegrateWW2(init=1,μ=1,ϵ=0.1,L=20,N=2^16,T=1,dt = 0.01,dealias=1,δ=0.01,m=-1/4,Ns=Ns)
		plt0=plot(problem0,var=[:surface,:fourier],T=0.6;compression=compression,label="")
		plt1=plot(problem1,var=[:surface,:fourier],T=0.6;compression=compression,label="")
		if !isnothing(name)
			savefig(plt0,"$name-N18.pdf");savefig(plt0,"$name-N18.svg");
			savefig(plt1,"$name-N16.pdf");savefig(plt1,"$name-N16.svg");
			if anim
				anim = @animate for t in LinRange(0,problem0.times.tfin,101)
					plt = plot([problem0,problem1];var=[:surface,:fourier],T=t,compression=compression,label=["N=2^18" "N=2^16"])
					ylims!(plt[1],(-0.5,1.1))			
				end	
				gif(anim, "$name.gif", fps = 15)
			end
		end
		return problem0,problem1,plt0,plt1

	elseif scenario in [7,7.1,7.2,7.3]
		if anim Ns=nothing else Ns=10 end
		p=round((scenario-7)*10);if p==0 p=2 end
		problem0,=IntegrateWW2(init=1,p=p,μ=1,ϵ=0.1,L=20,N=2^14,T=10,dt = 0.01,dealias=1,δ=0.01,m=-1,Ns=Ns)
		problem1,=IntegrateWW2(init=1,p=p,μ=1,ϵ=0.1,L=20,N=2^14,T=10,dt = 0.01,dealias=1,δ=0.002,m=-1,Ns=Ns)
		#problem2,=IntegrateWW2(init=1,p=p,μ=1,ϵ=0.1,L=20,N=2^14,T=2,dt = 0.01,dealias=1,δ=0.001,m=-1)
		plt0=plot(problem0;var=[:surface,:fourier],compression=compression,label="t=10")
		plot!(plt0,problem0,T=2;var=[:surface,:fourier],compression=compression,label="t=2")
		title!(plt0[1,1],"surface deformation")
		plt1=plot(problem1;var=[:surface,:fourier],compression=compression,label="t=10")
		plot!(plt1,problem1,T=2;var=[:surface,:fourier],compression=compression,label="t=2")
		title!(plt1[1,1],"surface deformation")
		if !isnothing(name)
			savefig(plt0,"$name-d01.pdf");savefig(plt0,"$name-d01.svg");
			savefig(plt1,"$name-d002.pdf");savefig(plt1,"$name-d002.svg");
			if anim
				anim = @animate for t in LinRange(0,problem0.times.tfin,101)
					plt = plot([problem0,problem1];var=[:surface,:fourier],T=t,compression=compression,label=["δ=0.01" "δ=0.002"])
					ylims!(plt[1],(-0.5,1.1))			
				end	
				gif(anim, "$name.gif", fps = 15)
			end
		end
		return problem0,problem1,plt0,plt1

	elseif scenario == 8
		blowups0=Real[];blowups1=Real[];
		K=100:20:800;iter=0;
		for k in K
			iter+=1;@info string("K=",k," (iteration ",iter,"/",length(K),")\n")
			~,blowup0=IntegrateWW2(init=2,K=k,μ=1,ϵ=0.15,L=20,N=2^14,T=10,dt = 0.001,dealias=1,δ=0,m=-1)
			~,blowup1=IntegrateWW2(init=2,K=k,μ=1,ϵ=0.2,L=20,N=2^14,T=10,dt = 0.001,dealias=1,δ=0,m=-1)
			push!(blowups0,blowup0);push!(blowups1,blowup1);
		end
		plt=scatter(K,[blowups0 blowups1],
				label=["ϵ=0.15" "ϵ=0.2"],
				marker=[:c :d],
				xlabel="K (in log scale)",
				ylabel="t (in log scale)",
				xscale=:log10,yscale=:log10)
		if !isnothing(name)
			savefig(plt,"$name.pdf");savefig(plt,"$name.svg");
			#@save(name,K,blowups0,blowups1)
		end
		return K,blowups0,blowups1,plt

	elseif scenario == 8.5
		K=100:20:800;
		blowups0=Real[];blowups1=Real[];iter=0;
		for k in K
			iter+=1;@info string("K=",k," (iteration ",iter,"/",length(K),")\n")
			N=2^14;d=N*π/k/20 # useful to define truncation at frequency K
			~,blowup0=IntegrateWW2(init=2,K=k,μ=1,ϵ=0.2,L=20,N=N,T=10,dt = 0.001,dealias=1,δ=0,m=-1)
			~,blowup1=IntegrateWW2(init=2,K=k,μ=1,ϵ=0.2,L=20,N=N,T=10,dt = 0.001,dealias=3/4*d-2,δ=0,m=-1)
			push!(blowups0,blowup0);push!(blowups1,blowup1);
		end
		plt=scatter(K,[blowups0 blowups1],
				label=["usual dealiasing" "adapted dealiasing"],
				marker=[:c :d],
				xlabel="K (in log scale)",
				ylabel="t (in log scale)",
				xscale=:log10,yscale=:log10)
		if !isnothing(name)
			savefig(plt,"$name.pdf");savefig(plt,"$name.svg");
			#@save(name,K,blowups0,blowups1)
		end
		return K,blowups0,blowups1,plt

	elseif scenario == 9
		Eps=10 .^(-1:0.02:0)
		blowups1=Real[];blowups2=Real[];blowups3=Real[];blowups4=Real[];iter=0;
		for eps in Eps
			iter+=1;@info string("ϵ=",eps," (iteration ",iter,"/",length(Eps),")\n")
			~,blowup1=IntegrateWW2(init=2,K=100,μ=1,ϵ=eps,L=20,N=2^14,T=10,dt = 0.001,dealias=1,δ=0,m=-1)
			~,blowup2=IntegrateWW2(init=2,K=200,μ=1,ϵ=eps,L=20,N=2^14,T=10,dt = 0.001,dealias=1,δ=0,m=-1)
			~,blowup3=IntegrateWW2(init=2,K=400,μ=1,ϵ=eps,L=20,N=2^14,T=10,dt = 0.001,dealias=1,δ=0,m=-1)
			~,blowup4=IntegrateWW2(init=2,K=800,μ=1,ϵ=eps,L=20,N=2^14,T=10,dt = 0.001,dealias=1,δ=0,m=-1)
			push!(blowups1,blowup1);push!(blowups2,blowup2);
			push!(blowups3,blowup3);push!(blowups4,blowup4);
		end
		plt=scatter(Eps,[blowups1 blowups2 blowups3 blowups4],
				label=["K=100" "K=200" "K=400" "K=800"],
				marker=[:c :d :ut :s],
				xlabel="ϵ (in log scale)",
				ylabel="t (in log scale)",
				xscale=:log10,yscale=:log10)
		if !isnothing(name)
			savefig(plt,"$name.pdf");savefig(plt,"$name.svg");
			#@save(name,Eps,blowups1,blowups2,blowups3,blowups4)
 		end
		return Eps,blowups1,blowups2,blowups3,blowups4,plt

	elseif scenario == 10
		K=400;N=2^14
		d=N*π/K/20 # useful to define truncation at frequency K
		problem0,=IntegrateWW2(init=2,K=K,μ=1,ϵ=0.2,L=20,N=N,T=2,dt = 0.001,dealias=1,δ=0,m=-1)
		problem1,=IntegrateWW2(init=2,K=K,μ=1,ϵ=0.2,L=20,N=N,T=2,dt = 0.001,dealias=3/4*d-2,δ=0,m=-1)

		plt0=plot(problem0, T=0.1;var=[:surface,:fourier],compression=compression,label="t=0.1")
		plot!(plt0,problem0,T=0.3;var=[:surface,:fourier],compression=compression,label="t=0.3")
		plot!(plt0,problem0,T=0.5;var=[:surface,:fourier],compression=compression,label="t=0.5")
		plot!(plt0[1,1],title="surface deformation")

		plt1=plot(problem1, T=0.1;var=[:surface,:fourier],compression=compression,label="t=0.1")
		plot!(plt1,problem1,T=0.3;var=[:surface,:fourier],compression=compression,label="t=0.3")
		plot!(plt1,problem1,T=0.5;var=[:surface,:fourier],compression=compression,label="t=0.5")
		plot!(plt1[1,1],title="surface deformation")
		if !isnothing(name)
			savefig(plt0,"$name-a.pdf");savefig(plt0,"$name-a.svg");
			savefig(plt1,"$name-b.pdf");savefig(plt1,"$name-b.svg");
			if anim
				anim = @animate for t in LinRange(0,problem0.times.tfin,101)
					plt = plot([problem0,problem1];var=[:surface,:fourier],T=t,compression=compression,label=["dealiasing" "low-pass filter"])
					ylims!(plt[1],(-0.5,1.1))			
				end	
				gif(anim, "$name.gif", fps = 15)

			end
		end
		return problem0,problem1,plt0,plt1

	elseif scenario == 11
		Eps=[0.01 0.0125 0.015 0.02 0.025 0.03 0.04 0.05 0.06 0.08 0.1 0.125 0.15 0.2 0.25 0.3 0.4 0.5 0.6 0.8]
		δc2=Real[];δc10=Real[];i=0;
		for eps in Eps
			i+=1
			#T=2
			δc_max=eps^2/2;δc_min=eps/500;
			for n=1:6
				@info string("iteration ",(i-1)*10+n,"/",10*length(Eps),"\n")
				@info string("critical ratio interval: ",(δc_min,δc_max))
				δc=sqrt(δc_max*δc_min)
				problem,blowuptime,blowup=IntegrateWW2(init=1,μ=1,ϵ=eps,L=20,N=2^20,T=2,dt = 0.005,dealias=1,δ=δc,m=-1,Ns=1)
				if blowup
					δc_min=δc
				else
					δc_max=δc
				end
			end
			push!(δc2,[δc_min,δc_max])
			#T=10
			δc_max=δc_min*1.5;
			for n=1:4
				@info string("iteration ",(i-1)*10+6+n,"/",10*length(Eps),"\n")
				@info string("critical ratio interval: ",(δc_min,δc_max))
				δc=sqrt(δc_max*δc_min)
				problem,blowuptime,blowup=IntegrateWW2(init=1,μ=1,ϵ=eps,L=20,N=2^20,T=10,dt = 0.005,dealias=1,δ=δc,m=-1,Ns=1)
				if blowup
					δc_min=δc
				else
					δc_max=δc
				end
			end
			push!(δc10,[δc_min,δc_max])
			#@save(name,Eps,δc2,δc10)
		end

		plt=plot()
		scatter!(plt,log10.(Eps[:]),[sqrt(δc10[i][2]*δc10[i][1]) for i in 1:length(Eps)],label="T=10",color=:1,marker=:c)
		scatter!(plt,log10.(Eps[:]),[sqrt(δc2[i][2]*δc2[i][1]) for i in 1:length(Eps)],label="T=2",color=:2,marker=:d)
		y2=OHLC[(δc2[i][1],δc2[i][1],δc2[i][2],δc2[i][2]) for i=1:length(Eps)]
		y10=OHLC[(δc10[i][1],δc10[i][1],δc10[i][2],δc10[i][2]) for i=1:length(Eps)]
		ohlc!(plt,log10.(Eps[:]),y10,label="",linecolor=:blue)
		ohlc!(plt,log10.(Eps[:]),y2,label="",linecolor=:orange)
		plot!(plt,log10.(Eps[:]),Eps[:].^2,label="δ=ϵ²",linecolor=:3)
		plot!(yscale=:log10,legend=:topleft,ylabel="δ (in log scale)",xlabel="ϵ (in log scale)")
		if !isnothing(name)
			savefig(plt,"$name.pdf");savefig(plt,"$name.svg");
			#@save(name,Eps,δc2,δc10)
		end
		return Eps,δc2,δc10,plt

	elseif scenario in [12,12.1,12.2,12.3]
		p=round((scenario-12)*10);if p==0 p=2 end

		Eps=[0.2,0.1,0.05]
		Delta=10 .^ (-2:0.1:0)
		Norms=Vector{Real}[];j=0;
		for eps in Eps
			norms=Real[];i=0;j+=1;
			for delta in Delta
				i+=1
				@info string("ϵ=",eps,", δ=",delta," (iteration ",(j-1)*length(Delta)+i,"/",length(Eps)*length(Delta),")")
				problem,=IntegrateWW2(init=1,p=p,μ=1,ϵ=eps,L=20,N=2^12,T=10,dt = 0.01,dealias=1,δ=delta,m=-1)
				param = (μ=1,ϵ=eps,L=20,N=2^12,T=10,dt = 0.01)
				model = WaterWaves(param; dealias = 1,method=1,maxiter=20)
				problem0 = Problem(model, problem.initial, param)
				solve!( problem0 )
				(η0,v0,x0)=solution(problem0)
				(η,v)=solution(problem,x=x0)
				push!(norms,norm(η-η0,Inf))
			end
			push!(Norms,norms)
		end

		plt=plot()
		markers=[:c :d :ut :s :r :lt :pent :hex :hep :oct :+ :x]
		for i in 1:length(Eps)
			scatter!(plt,Delta,Norms[i],label="ϵ=$Eps[i]",marker=markers[i])
		end
		plot!(plt,xlabel="δ",ylabel="error",legend=:topleft)
		if !isnothing(name)
			savefig(plt,"$name.pdf");savefig(plt,"$name.svg");
			plot!(plt,xscale=:log10,yscale=:log10)
			plot!(plt,xlabel="δ (in log scale)",ylabel="error (in log scale)")
			savefig(plt,"$name-log.pdf");savefig(plt,"$name-log.svg");
			#@save(name,Eps,Delta,Norms)
		end
		return Eps,Delta,Norms,plt

	else
		error("The first argument must be as specified in the documentation (basically between 1 and 12)")
	end

end
nothing
