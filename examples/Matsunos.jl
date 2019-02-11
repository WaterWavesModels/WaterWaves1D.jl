# # Comparison of Matsuno deep water problems
#using DeepWaterModels
include("../src/dependencies.jl")

#----
param = ( ϵ  = 1/2,
          N  = 2^12,
          L  = 10,
          T  = 5,
          dt = 0.001)

init     = BellCurve(param,2.5)

model0   = Matsuno_mod_naive(param)
problem0 = Problem(model0, init, param);

model1   = Matsuno_naive(param)
problem1 = Problem(model1, init, param);

model2   = Matsuno(param)
problem2 = Problem(model2, init, param);

using DelimitedFiles

open("reference.dat", "w") do io
    writedlm(io, problem1.data.U[end])
end

problems = [ problem0, problem1, problem2 ]


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
