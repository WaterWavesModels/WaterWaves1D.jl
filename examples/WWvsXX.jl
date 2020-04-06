# # Two water problems
#
#md # [`notebook`](@__NBVIEWER_ROOT_URL__notebooks/two_problems.ipynb)
#
#using ShallowWaterModels
include("../src/dependencies.jl")

#----
param = ( μ  = .2,
			ϵ  = .5,
        	N  = 2^10,
            L  = 20,
            T  = 5,
            dt = 0.01,
			θ = 1,
			α = 1,
			a = -1/3, b = 1/3)
g(x) = exp.(-x.^2)
z(x) = 0*exp.(-x.^2)

init     = Init(g,z)

model0  = WhithamGreenNaghdiSym(param;SGN=true,iterate=true,precond=true)
solver0   = RK4(param,model0)
problem0 = Problem(model0, init, param; solver = solver0);

model1    = WhithamGreenNaghdi(param;SGN=true,iterate=true)
solver1   = RK4(param,model1)
problem1 = Problem(model1, init, param; solver = solver1);

# model2    = WhithamGreenNaghdi(param;SGN=true,iterate=true,precond=false)
# solver2   = RK4(param,model2)
# problem2 = Problem(model2, init, param; solver = solver2);



model2    = WaterWaves(param)
solver2   = RK4(param,model2)
problem2 = Problem(model2, init, param; solver = solver2);

#----


problems = [ problem0, problem1, problem2 ]

for problem in problems
	print("\nNow solving the model ",problem.model.label,"\n")
   	@time solve!( problem )
end

p = plot(layout=(2,1))

for problem in problems
   	fig_problem!( p, problem)
end
display(p)
savefig("WWvsXX.png"); nothing # hide

#create_animations(problems)

#----
#md # ![](two_problems.png)
