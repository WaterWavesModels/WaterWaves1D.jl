# # Two water problems
#
#md # [`notebook`](@__NBVIEWER_ROOT_URL__notebooks/two_problems.ipynb)
#
#using DeepWaterModels
include("../src/dependencies.jl")

#----
param = ( μ  = 1/2,
			ϵ  = 0.1,
        	N  = 2^11,
            L  = 20,
            T  = 15,
            dt = 0.001,
			theta = 2)

init     = BellCurve(param)

model0    = WaterWaves(param)
solver0   = RK4(param,model0)
problem0 = Problem(model0, init, param, solver0);

model1    = fdBoussinesq_1(param)
solver1   = RK4(param,model1)
problem1 = Problem(model1, init, param, solver1);

model2  = Boussinesq(param)
solver2   = RK4(param,model2)
problem2 = Problem(model2, init, param, solver2);


#----

p = plot(layout=(2,1))

problems = [ problem0, problem1, problem2 ]

for problem in problems
	print("\nNow solving the model ",problem.model.label,"\n")
   	@time solve!( problem )
   	fig_problem!( p, problem )

end

#savefig("two_problems.png"); nothing # hide
display(p)
#----
#md # ![](two_problems.png)
