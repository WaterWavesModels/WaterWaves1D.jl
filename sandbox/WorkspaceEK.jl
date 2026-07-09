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
include("../src/models/EulerKorteweg.jl")
using FFTW
#---- parameters
param = (
    κ = 0.001,
    ϵ = 1,
    N = 2^12,     # number of collocation points
    L = π,    # size of the mesh (-L,L)
    T = 0.5,        # final time of computation
    dt = 0.0001, # timestep
    ns = 50,       # data stored every ns step
);
#---- initial data
η(x) = -exp.(-abs.(x) .^ 3);
u(x) = 0.1 * (x .+ 1) .* exp.(-x .^ 2);
init = Init(η, u);

#---- building problems
model1 = EulerKorteweg(param)
model2 = EulerKorteweg_GP(param)
model3 = EulerKorteweg_Grenier(param)
model = SaintVenant(param; dealias = 1)
model0 = SaintVenant(param; dealias = 1, smooth = 0.1)


problems = Problem[]
#push!(problems, Problem(model2, init, param, solver= Euler(model2),label="GPEuler" ))
#push!(problems, Problem(model2, init, param, solver= EulerExp(model2),label="GPEulerSym" ))
push!(problems, Problem(model3, init, param, label = "Grenier"))
#push!(problems, Problem(model1, init, param, label = "EK" ))
push!(problems, Problem(model2, init, param, label = "GP"))
push!(problems, Problem(model, init, param, label = "SV"))
push!(problems, Problem(model0, init, param, label = "SVsmooth"))

#push!(problems, Problem(model, init, param, solver=solver2, label="exponential naive" ))
#push!(problems, Problem(model, init, param, solver=solver3, label="exponential" ))


#---- computation
for problem in problems
    @time solve!(problem)
end
#solve!(problems[3])
#---- visualization
#include("../src/Figures.jl")
#using Plots
#gr()
using Plots
plot(problems, T = 0.3, var = [:surface])

(η, v, x, t) = solution(problems[1])
(η2, v2, x2, t) = solution(problems[2])
(η3, v3, x3, t) = solution(problems[3])
(η4, v4, x4, t) = solution(problems[4])

using Plots;gr()
plot(x, [η - η2])
plot(x, [η - η3 η - η4])
plot(x, [η4 - η3])
