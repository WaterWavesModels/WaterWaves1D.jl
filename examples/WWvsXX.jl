# # Two water problems
#
#md # [`notebook`](@__NBVIEWER_ROOT_URL__notebooks/two_problems.ipynb)
#
#using ShallowWaterModels
include("../src/dependencies.jl")

#----
param = ( μ  = 0.1,
			ϵ  = 0.1,
        	N  = 2^11,
            L  = 20,
            T  = 15,
            dt = 0.001,
			θ = 1,
			α = 1,
			a = -1/3, b = 1/3)

init     = Bellcurve(param)

model0    = WaterWaves(param)
solver0   = RK4(param,model0)
problem0 = Problem(model0, init, param, solver0);

model1    = Boussinesq(param)
solver1   = RK4(param,model1)
problem1 = Problem(model1, init, param, solver1);

model2  = fdBoussinesq(param)
solver2   = RK4(param,model2)
problem2 = Problem(model2, init, param, solver2);


#----


problems = [ problem0, problem1, problem2 ]

for problem in problems
	print("\nNow solving the model ",problem.model.label,"\n")
   	@time solve!( problem )
end

p = plot(layout=(2,1))

for problem in problems
   	fig_problem!( p, problem )
end
savefig("WWvsXX.png"); nothing # hide
display(p)
#create_animations(problems)

#----
#md # ![](two_problems.png)
