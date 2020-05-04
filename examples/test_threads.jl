# # Comparing several models for water waves
#
#md # [`notebook`](@__NBVIEWER_ROOT_URL__notebooks/WWvsXX.ipynb)
#
using ShallowWaterModels;
using Test
# if there is an error at this step, try commenting the lines above and commenting out the line below
#include("../src/dependencies.jl")

#---- parameters
param = ( μ  = .1,
			ϵ  = 1,
        	N  = 2^11, 	# number of collocation points
            L  = 10,	# size of the mesh (-L,L)
            T  = 5,		# final time of computation
            dt = 0.001, # timestep
				);
#---- initial data
g(x) = exp.(-abs.(x).^4);
z(x) = 0*exp.(-x.^2);
init = Init(g,z);

#---- models to compare
models=[]
push!(models,WaterWaves(param))
push!(models,PseudoSpectral(param;order=2,dealias=1,lowpass=1/100))
push!(models,PseudoSpectral(param;order=3,dealias=1,lowpass=1/100))

#---- Initialize
function Problems(models)
	problems = []
	for model in models
		push!(problems, Problem(model, init, param) )
	end
	return problems
end

#---- visualization
function Plot_problems(problems)
	p=plot(layout=(2,1))
	for problem in problems
		fig_problem!( p, problem)
	end
	display(p)
end

#---- solving

# using Distributed
# #addprocs(3)
#
# problems=Problems(models)
# @time @async for problem in problems
# 	solve!( problem )
# end

#---- Method 1
@info "Method 1 with threads"
problems=Problems(models)
@time solve!(problems)
p1=(last(problems[1].data.U), last(problems[2].data.U) ,last(problems[3].data.U))

#---- Method 2
@info "Method 2 with threads"
problems=Problems(models)
using Base.Threads
@time @threads for problem in problems
	solve!( problem )
end
p2=(last(problems[1].data.U), last(problems[2].data.U) ,last(problems[3].data.U))

#---- Method 3
@info "Consecutive solving"
problems=Problems(models)
@time for problem in problems
	solve!( problem )
end
p3=(last(problems[1].data.U), last(problems[2].data.U) ,last(problems[3].data.U))

#---- Tests
@testset "Final" begin
	@test p1==p2
	@test p2==p3
end
