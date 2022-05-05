using Plots
using WaterWaves1D

param = (
    μ  = 1, 
    ϵ  = 1/4, 
    N  = 2^11,
    L  = 10, 
    T  = 5, 
    dt = 0.01)

z(x) = exp.(-abs.(x).^4)
v(x) = zero(x)
init = Init(z,v)

model1 = WaterWaves(param, tol = 1e-15)
model2 = WWn(param;n=2,dealias=1,δ=1/10)

problem1 = Problem(model1, init, param)
problem2 = Problem(model2, init, param)

solve!([problem1 problem2]; verbose=false)

#p = plot(layout=(1,2))
#plotdifferences!(p[1,1], problem1, problem2)
plotdifferences([problem1, problem2])
