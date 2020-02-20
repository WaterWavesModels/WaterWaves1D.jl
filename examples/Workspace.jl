# # Workspace
#
#md # [`notebook`](@__NBVIEWER_ROOT_URL__notebooks/two_problems.ipynb)
#
#using ShallowWaterModels
include("../src/dependencies.jl")

#----
param = ( μ  = 1/2,
			ϵ  = 1,
        	N  = 2^14,
            L  = 5,
            T  = 1.5,
            dt = 0.001)

#init     = BellCurve(merge(param,(θ = 1,)))
#init = Init((η = x-> 2 .^(-abs.(x.^2)) , v = x->0 .*x))
η = x-> 10*x.*2 .^(-abs.(x.^4))
v = x->-η(x)
init = Init(η,v)


models = []
push!(models,fdBoussinesq(merge(param,(α=1,))))
push!(models,fdBoussinesq(merge(param,(α=1/2,))))
push!(models,fdBoussinesq(merge(param,(α=2/5,))))

problems = []
for model in models
	push!(problems, Problem(model, init, param))
end


#----

for problem in problems
	print("\nNow solving the model ",problem.model.label,"\n")
   	@time solve!( problem )
end

p = plot(layout=(2,1))

for problem in problems
   	fig_problem!( p, problem, 0.2)
end
#plot!(xlims=[(-1,-0.5) (-3000,3000 )])
display(p)

savefig("final.pdf"); nothing # hide
#create_animations(problems)
