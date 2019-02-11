include("../src/dependencies.jl")

param = ( ϵ  = 1/2,
          N  = 2^12,
          L  = 10,
          T  = 5,
          dt = 0.001 )

init    = BellCurve(param,2.5)
model   = Matsuno(param)
problem = Problem(model, init, param);

p = plot(layout=(2,1))

@time solve!( problem )
fig_problem!( p, problem )
savefig("test.png")
