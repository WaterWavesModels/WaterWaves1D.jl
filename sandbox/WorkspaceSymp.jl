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
#include("../src/solvers/StoermerVerlet.jl")


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
WW_model = WaterWaves(param)
WW_model = WWn(param;n=2,dealias=1,δ=0.1)

RK4_problem = Problem(WW_model, init, param, solver = RK4(WW_model), label = "Runge-Kutta 4")
Euler_problem = Problem(WW_model, init, param, solver = Euler(WW_model), label = "Explicit Euler")
Symplectic_problem = Problem(WW_model, init, param, solver = EulerSymp(WW_model, Niter = 5, implicit = 2), label = "Symplectic Euler")
StoermerVerlet_problem = Problem(WW_model, init, param, solver = StoermerVerlet(WW_model, Niter = 5, implicit = 2), label = "Stoermer-Verlet")

#Symplectic_problem = Problem(WW_model, init, param, solver = EulerSymp(WW_model,Niter=5,implicit=2))

for problem in [RK4_problem, Euler_problem, Symplectic_problem, StoermerVerlet_problem]
    solve!(problem)
end

using Plots
plot([RK4_problem, Euler_problem, Symplectic_problem, StoermerVerlet_problem])

plot((RK4_problem, Euler_problem))
plot((RK4_problem, Symplectic_problem))
plot((RK4_problem, StoermerVerlet_problem))

energy_diff(RK4_problem)
energy_diff(Euler_problem)
energy_diff(Symplectic_problem)
energy_diff(StoermerVerlet_problem)


#solve!(problems[3])
#---- visualization
#include("../src/Figures.jl")
#using Plots
#gr()
(η, v, x, t) = solution(problems[1])
(η2, v2, x2, t) = solution(problems[2])
(η3, v3, x3, t) = solution(problems[3])
(η4, v4, x4, t) = solution(problems[4])

using Plots;gr()
plot(x1, [η - η2])
plot(x1, [η - η2 η - η3])
