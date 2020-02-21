# # Workspace
#
#md # [`notebook`](@__NBVIEWER_ROOT_URL__notebooks/two_problems.ipynb)
#
#using ShallowWaterModels
include("../src/dependencies.jl")

#----
param = ( μ  = 1/2,
			ϵ  = 1,
        	N  = 2^17,
            L  = 3,
            T  = 1,
            dt = 0.0001)

#init     = BellCurve(merge(param,(θ = 1,)))
#init = Init((η = x-> 2 .^(-abs.(x.^2)) , v = x->0 .*x))
η = x -> 2 .^(-abs.((x).^4))
v = x -> -10*x.*η(x)
init = Init(η,v)


models = []
#push!(models,fdBoussinesq(merge(param,(α=1,))))
push!(models,fdBoussinesq(merge(param,(α=1/2,))))
#push!(models,fdBoussinesq(merge(param,(α=2/5,))))

problems = []
for model in models
	push!(problems, Problem(model, init, param))
end


#----

for problem in problems
	print("\nNow solving the model ",problem.model.label,"\n")
   	@time solve!( problem )
end

#----
p = plot(layout=(2,1))
for problem in problems
   	fig_problem!( p, problem, .92)
end
plot!(xlims=[(0.275,0.285) (-100000,100000 )])
display(p)

savefig("final.pdf"); nothing # hide
#create_animations(problems)
#p=problems[1]
#plot(p.mesh.x,mapfro(p.model,p.data.U[100]))
