# # Workspace
#
#md # [`notebook`](@__NBVIEWER_ROOT_URL__notebooks/two_problems.ipynb)
#
using WaterWaves1D, Plots, FFTW, Statistics;
gr();
#include("../src/dependencies.jl")

include("../src/models/WaterWaves.jl")
include("../src/models/DeepWaterWaves.jl")
include("../src/models/WWn.jl")

#include("../src/models/SerreGreenNaghdi.jl")
#include("../src/models/Nonhydrostatic.jl")
#include("../src/models/SquareRootDepth.jl")

include("../src/Figures.jl")


# if there is an error at this step, try commenting the lines above and commenting out the line below
#include("../src/dependencies.jl")

#---- parameters
param = (
    μ = 10^2,
    ϵ = 0.1/10^1,
    N = 2^8, # number of collocation points
    L = 5,# size of the mesh (-L,L)
    T = 1,# final time of computation
    dt = 0.001,  # timestep
);
param2 = (
    μ = Inf,
    ϵ = 0.1,
    N = 2^8, # number of collocation points
    L = 5,# size of the mesh (-L,L)
    T = 1,# final time of computation
    dt = 0.001,  # timestep
);

#---- initial data
ζ(x) = exp.(-abs.(x) .^ 4);
v2(x) = zero(x);
init = Init(ζ, v2);
# ζ et v sont des fonctions.
# Init les mets sous une forme utilisable par le programme


#---- models to compare
models = []
push!(models, WaterWaves(param; dealias = 0))  # dealias = 1 pour dealiasing
#push!(models, WaterWaves(param; dealias = 0))  # dealias = 1 pour dealiasing
#push!(models, WWn(param; n=2,dealias = 1))  # dealias = 1 pour dealiasing
push!(models, WaterWaves(param2; dealias = 0))  # dealias = 1 pour dealiasing
#push!(models, WaterWaves(param; dealias = 0))  # dealias = 1 pour dealiasing
#push!(models, WWn(param2; n=2,dealias = 1))  # dealias = 1 pour dealiasing


problems = []
for model in models
    push!(problems, Problem(model, init, param))
end

#---- computation
for problem in problems
    @time solve!(problem)
end

#---- visualization
plt = plot_solution(problems; fourier = false)
plot_difference(problems; fourier = false)
#savefig(plt,"plot.pdf") # ou "plot.png" par exemple

# pour visualiser autre chose, les différences par exemple
#(η1,v1,x1,t)=solution(problems[1])
#(η2,v2,x2,t)=solution(problems[2],x=x1)
#(η3,v3,x3,t)=solution(problems[3],x=x1)

#plt = plot(x1,η1-η2,label="diff 1")
#plot!(plt,x1,η1-η3,label="diff 2")
#plot!(xlims=(0,10))
#plot!(ylims=(-0.1,0.1))

#plot_solution(problems;fourier=false)
# pour faire un film
#anim = create_animation(problems;fourier=false,ylims=(-0.5,1.1))
#gif(anim, "SGNvsWGN.gif", fps=15)
#mov(anim, "SGNvsWGN.mov", fps=15)
#mp4(anim, "SGNvsWGN.mp4", fps=15)
