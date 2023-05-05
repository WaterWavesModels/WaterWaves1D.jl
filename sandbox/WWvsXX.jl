# # Comparingseveral models for water waves
#
#md # [`notebook`](@__NBVIEWER_ROOT_URL__notebooks/WWvsXX.ipynb)
#
using WaterWaves1D;
#include("../src/models/WaterWaves.jl")
#include("../src/models/PseudoSpectral.jl")

# if there is an error at this step, try commenting the lines above and commenting out the line below
#include("../src/dependencies.jl")

#---- parameters
param = ( μ  = 0.1,
			ϵ  = 0.05,
        	N  = 2^12, 	# number of collocation points
            L  = 15,	# size of the mesh (-L,L)
            T  = 10,		# final time of computation
            dt = 0.001, # timestep
			ns = 50,   	# data stored every ns step
				);
#---- initial data
g(x) = exp.(-abs.(x).^4);
z(x) = 0*exp.(-x.^2);
init = Init(g,z);

#---- models to compare
models=AbstractModel[]
push!(models,WaterWaves(param))
push!(models,Whitham(param))
push!(models,Whitham(param;BBM=true))
push!(models,WhithamBoussinesq(param))


#push!(models,PseudoSpectral(param;order=2,dealias=1,lowpass=1/100))
#push!(models,PseudoSpectral(param;order=3,dealias=1,lowpass=1/100))

problems = Problem[]
for model in models
	push!(problems, Problem(model, init, param) )
end

#---- computation
for problem in problems
	solve!( problem )
end

#---- visualization
#include("../src/Figures.jl")
using Plots
gr()

plt = plot(problems;xlims=(0,15),legend=:topleft)
display(plt)
savefig("WWvsXX.pdf")

plot(problems[1:2];color=:2,xlims=(0,15),var=:difference)
plot!(problems[[1 3]];color=:3,xlims=(0,15),var=:difference)
plot!(problems[[1 4]];color=:4,xlims=(0,15),var=:difference)
savefig("differences.pdf")


anim = @animate for time in LinRange(0,param.T,101)
    plot(problems,T=time,ylims=(-0.2,1))
end
gif(anim,"WWvsXX.gif",fps=15)