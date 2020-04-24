# # Two deep water problems
#
#md # [`notebook`](@__NBVIEWER_ROOT_URL__notebooks/two_problems.ipynb)
#
#using ShallowWaterModels
include("../src/dependencies.jl")

#----
param = ( ϵ  = 1/2,
        	N  = 2^12,
            L  = 10,
            T  = 5,
            dt = 0.001,
			θ = 2,
			s = 2,
			freq =1)

init     = Bellcurve(param)

model1    = CGBSW(param)
solver1   = RK4(param,model1)
problem1 = Problem(model1, init, param; solver = solver1);

model2  = Matsuno(param)
solver2   = RK4(param,model2)
problem2 = Problem(model2, init, param; solver = solver2);

problems = [ problem1, problem2 ]

#----

p = plot(layout=(2,1))


for problem in problems
	print("\nNow solving the model ",problem.model.label,"\n")
   	@time solve!( problem )
   	fig_problem!( p, problem )

end

#savefig("two_problems.png"); nothing # hide
display(p)
#----
#md # ![](two_problems.png)
