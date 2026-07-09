# # Workspace
#
#md # [`notebook`](@__NBVIEWER_ROOT_URL__notebooks/two_problems.ipynb)
#
#using WaterWaves1D,ProgressMeter,Plots,FFTW,Statistics;gr()
#include("../src/dependencies.jl")
#
#include("../src/models/PseudoSpectral.jl")
#include("../src/models/Spectral.jl")
#include("../src/models/WaterWaves.jl")
#
#include("../src/solvers/Euler.jl")
#include("../src/solvers/EulerSymp.jl")


# if there is an error at this step, try commenting the lines above and commenting out the line below
#include("../src/dependencies.jl")

#---- parameters
param = (
    μ = 0.25,
    ϵ = 0.1,
    N = 2^8,     # number of collocation points
    L = π,    # size of the mesh (-L,L)
    T = 50,        # final time of computation
    dt = 0.001, # timestep
    ns = 50,       # data stored every ns step
);
#---- initial data
g(x) = exp.(-abs.(x) .^ 4);
z(x) = 0 .* (x .+ 1) .* exp.(-x .^ 2);
init = Init(g, z);

#---- building problems
model = BBM(param)
solver1 = Euler_naive()
solver2 = EulerExp_naive()
solver3 = EulerExp(model)


problems = Problem[]
push!(problems, Problem(model, init, param, label = "RK4"))
push!(problems, Problem(model, init, param, solver = solver1, label = "Euler"))
push!(problems, Problem(model, init, param, solver = solver2, label = "exponential naive"))
push!(problems, Problem(model, init, param, solver = solver3, label = "exponential"))


#---- computation
for problem in problems
    @time solve!(problem)
end
#solve!(problems[3])
#---- visualization
#include("../src/Figures.jl")
#using Plots
#gr()

plot(problems, T = 50)

(η, v, x, t) = solution(problems[1])
(η2, v2, x2, t) = solution(problems[2])
(η3, v3, x3, t) = solution(problems[3])
(η4, v4, x4, t) = solution(problems[4])

using Plots;gr()
plot(x, [η - η2])
plot(x, [η - η3 η - η4])
plot(x, [η4 - η3])
