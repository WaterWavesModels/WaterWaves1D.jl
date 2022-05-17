# #
# Reproduce the figures in the Appendix of
# [Duchêne](https://www.ams.org/open-math-notes/omn-view-listing?listingId=111309)
# comparing solutions to the water waves system with predictions of
# the Serre-Green-Naghdi and Whitham-Green-Naghdi systems
# #
export Integrate,Figures
using WaterWaves1D,Plots,LinearAlgebra;
#using JLD #(uncomment if using @save)

#----
"""
	Integrate(scenario;kwargs)

Integrate in time the water waves system as well as
the Serre-Green-Naghdi, Whitham-Green-Naghdi and Isobe-Kakinuma models
with an initial datum depending on the argument (`scenario`).
If `scenario==1`, the initial datum is a heap of water with zero velocity.
If `scenario==2`, the initial datum is generates at first order a right-going wave.


Other arguments are optional:
- `μ` is the shallowness parameter,
- `ϵ` in the nonlinearity parameter,
- `p` is the power in the initial datum `η = exp(-|x|^p)`
- `N` the number of collocation points,
- `L` the half-length of the mesh,
- `T` the final time of integration,
- `dt` the timestep,
- `dealias`: dealiasing with Orlicz rule `1-dealias/(dealias+2)` (default is `0`, i.e. no dealiasing),
- `iterate`: use GMRES if `true`, and LU decomposition otherwise (default is `true`),
- `precond` gives some choice in the preconditioner for GMRES (default is `true` for the default preconditioner),
- `method` and `maxiter` are the corresponding options for generating the initial data in the water waves model; see doc of `WaterWaves`
- `name`: a string used to save the figures.

Return `(Diff,problems,plt)` where `Diff` is a measure of the errors, `problems` contains all the information and `plt` is a plot of the final time solution.
"""
function Integrate(scenario;μ=.1,ϵ=.1,p=2,N=2^10,L=15,T=10,dt=0.001,dealias=0,iterate=true,precond=true,method=1,maxiter=100,name=nothing)
	#---- preparation

	if !isnothing(name) ns=floor(Int,max(1,T/dt/100)) else ns=1 end

	if typeof(dealias)==Int64 dealias=dealias*ones(Int64,4) end

	# the set of parameters as indicated by the user
	param = ( μ  = μ, ϵ  = ϵ,
				N  = N, L  = L,
	            T  = T, dt = dt,
				ns=ns )

	# generate the initial data
	mesh=Mesh(param)
	k=mesh.k

	η = exp.(-abs.(mesh.x ).^p)
	if scenario == 1
		v = 0*η
	elseif scenario == 2
		v = 2/ϵ*(sqrt.(1 .+ ϵ*η) .- 1)
	else
		@error("the argument must be 1 or 2")
	end
	init     = Init(mesh,η,v)

	# construct the models to be solved
	modelWW	 = WaterWaves(param;dealias=dealias[1],method=method,maxiter=maxiter)
	modelSGN = SerreGreenNaghdi(param;ktol=0, gtol=1e-14, iterate=iterate, precond = precond, dealias = dealias[2])
	modelWGN = WhithamGreenNaghdi(param;SGN=false, ktol=0, gtol=1e-14, iterate=iterate, precond = precond, dealias = dealias[3])
	modelIK2 = IsobeKakinuma(param; ktol=0, gtol=1e-14, iterate=iterate, precond = precond, dealias = dealias[4])

	# construct the associated problems
	problems = []
	push!(problems, Problem(modelWW, init, param) )
	push!(problems, Problem(modelSGN, init, param) )
	push!(problems, Problem(modelWGN, init, param) )
	push!(problems, Problem(modelIK2, init, param) )



	#---- computation
	for problem in problems
		solve!( problem )
	end

	#---- analysis
	# Errors at final time
	(ηww,vww,xww,t)=solution(problems[1])
	(ηSGN,vSGN)=solution(problems[2],x=xww)
	(ηWGN,vWGN)=solution(problems[3],x=xww)
	(ηIK2,vIK2)=solution(problems[4],x=xww)

	Diff = [norm(ηww-ηSGN,Inf) norm(ηww-ηWGN,Inf) norm(ηww-ηIK2,Inf)]
	@info string("Error of SGN, WGN, IK2: ",round.(Diff, sigdigits=3))

	#---- plots
	# solutions at final time
	plt = plot(problems)
	display(plt)


	if !isnothing(name)
		#save(problem,name);
		#@save(name,xww,ηww,ηSGN,ηWGN,ηIK2,scenario,μ,ϵ,p,N,L,T,dt,dealias,iterate,precond,name)

		if scenario == 1
			xlims=(0,last(xww))
			ind=range(1,step=length(xww)÷32,length=32);
		else
			xlims=(first(xww),last(xww))
			ind=range(1,step=length(xww)÷16,length=16);
		end

		# Plot and save the solutions at final time
		plot(problems;fourier=false)
		plot!(xlims=xlims,legend=:topleft)
		savefig(string(name,".pdf"));

		# Plot and save the differences at final time
		plot(xww, [ηww-ηSGN ηww-ηWGN ηww-ηIK2],
		 	color=[:2 :3 :4],primary=false)
		scatter!(xww[ind], [ηww[ind]-ηSGN[ind] ηww[ind]-ηWGN[ind] ηww[ind]-ηIK2[ind]],
		 	color=[:2 :3 :4],marker=([:d :h :ut], 8, Plots.stroke(1, :gray)),#linewidth=0,
		 	xlims=xlims,
		   title=string("difference (surface deformation) at t=",t),
		   label=["Serre-Green-Naghdi" "Whitham-Green-Naghdi" "Isobe-Kakinuma"])
		savefig(string(name,"-differences.pdf"));

		# Generate and save an animation of the evolution in time
		anim = @animate for time in LinRange(0,param.T,100)
			plot(problems; T = time, xlims=xlims,ylims=(-0.5,1))
		end
		gif(anim,string(name,".gif"))

	end

	return Diff,problems,plt
end

#----
"""
	Figures()

Reproduce the figures in Section I.5 of the monograph
[Many Models for Water Waves](https://www.ams.org/open-math-notes/omn-view-listing?listingId=111309)
"""
function Figures()
	MU=[1 0.1 0.01]
	EPS=[0.5 0.25 0.125]

	DiffSGN=zeros(3,3);DiffWGN=zeros(3,3);DiffIK2=zeros(3,3);

	for i in 1:3
		for j in 1:3
			mu=MU[i];eps=EPS[j]
			if mu == 1. mu = 1 end
			name=string("mu",mu,"eps",eps)
			name=replace(name, '.' => "")
			@info string("Save plots and animation with the name ",name)
			Diff,=Integrate(1;p=3,μ=mu,ϵ=eps,N=2^10,T=10,L=15,dt=0.001,dealias=0,maxiter=100,name=name);
			DiffSGN[i,j]=Diff[1];DiffWGN[i,j]=Diff[2];DiffIK2[i,j]=Diff[3];
			#@save("MM4WW",DiffSGN,DiffWGN,DiffIK2)
		end
	end
end
nothing
