# # Two deep water problems
#
#md # [`notebook`](@__NBVIEWER_ROOT_URL__notebooks/two_problems.ipynb)
#
#using DeepWaterModels
include("../src/dependencies.jl")

#----
param = ( ϵ  = 1/2,
        	N  = 2^12,
            L  = 10,
            T  = 5,
            dt = 0.001)

init     = BellCurve(param)
solver   = RK4(param,2)

model1    = CGBSW(param)
problem1 = Problem(model1, init, param, solver);

model2  = Matsuno(param)
problem2 = Problem(model2, init, param, solver);


#----

p = plot(layout=(2,1))

problems = [ problem1, problem2 ]

for problem in problems
	print("\nNow solving the model ",problem.model.label,"\n")
   	@time solve!( problem )
   	fig_problem!( p, problem )

end

#savefig("two_problems.png"); nothing # hide
display(p)
#----
#md # ![](two_problems.png)
