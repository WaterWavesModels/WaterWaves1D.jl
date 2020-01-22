# # Two water problems
#
#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/notebooks/two_problems.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/notebooks/two_problems.ipynb)
#
using ShallowWaterModels

#----
param = ( μ  = 1/20,
          ϵ  = 1/2,
          N  = 2^11,
          L  = 20,
          T  = 15,
          dt = 0.001,
          theta = 2)

init = BellCurve(param)

model1 = fdBoussinesq_1b(param)
solver1 = RK4(param,model1)
problem1 = Problem(model1, init, param, solver1);

model2 = fdBoussinesq_1(param)
solver2 = RK4(param,model2)
problem2 = Problem(model2, init, param, solver2);

#----

p = plot(layout=(2,1))

problems = [ problem1, problem2 ]

for problem in problems
	print("\nNow solving the model ",problem.model.label,"\n")
   	@time solve!( problem )
   	fig_problem!( p, problem )

end

#nb # display(p)
#----
savefig("two_problems.png")
#md # ![](two_problems.png)
