# # Saved Experiments
#
#md # [`notebook`](@__NBVIEWER_ROOT_URL__notebooks/two_problems.ipynb)
#
#using ShallowWaterModels
include("../src/dependencies.jl")

#---- Wavebreaking for the quasilinear Boussinesq-Whitham equation
param = ( μ  = 1/2,
			ϵ  = 1,
        	N  = 2^18,
            L  = 3,
            T  = .04,
            dt = 0.00001)
η = x -> 2 .^(-abs.((x).^4))
v = x -> -100*x.*η(x)
init = Init(η,v)
model = fdBoussinesq(merge(param,(α=1/2,)))
problem = Problem(model, init, param)
@time solve!( problem )
p = plot(layout=(2,1))
fig_problem!( p, problem, .03898)
plot!(xlims=[(0.0018,0.0025) (-100000,100000 )])
display(p)


#---- Wavebreaking for the quasilinear Boussinesq-Whitham equation ?
param = ( μ  = 1/2,
			ϵ  = 1,
        	N  = 2^17, # you need to add more modes
            L  = 3,
            T  = 1,
            dt = 0.0001) # you need to lower this
η = x -> 2 .^(-abs.((x).^4))
v = x -> -10*x.*η(x)
init = Init(η,v)
model = fdBoussinesq(merge(param,(α=1/2,)))
problem = Problem(model, init, param)
@time solve!( problem )
p = plot(layout=(2,1))
fig_problem!( p, problem, .92)
plot!(xlims=[(0.275,0.285) (-100000,100000 )])
display(p)
