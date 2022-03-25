# Plot recipes

There some [recipes](https://docs.juliaplots.org/latest/recipes/) available to visualize
problems solution. You need to import the package [Plots.jl](https://github.com/JuliaPlots/Plots.jl).


```@example surface
using Plots
using WaterWaves1D

param = ( μ = 1, ϵ = 1/4, N = 2^10, L = 10, T = 5, dt = 0.01 )
z(x) = exp.(-abs.(x).^4)
v(x) = zero(x)
init = Init(z,v)

model1 = WaterWaves(param) # The water waves system
model2 = WWn(param;n=2,dealias=1,δ=1/10) # The quadratic model (WW2)

problem1 = Problem(model1, init, param) ;
problem2 = Problem(model2, init, param) ;

solve!([problem1 problem2]; verbose=false);

plot(problem1)
plot!(problem2; legend = :bottomright)
```

```@example surface
plot([problem1, problem2])
```

```@example fourier
using Plots
using WaterWaves1D

η(x) = exp.(-x.^2)
v(x) = zero(x)   
init = Init(η,v)

param1 = ( ϵ = 1/4, μ = Inf, N = 2^8, L = 2*π, T = 5, dt = 0.001 )

problem1 = Problem( WWn(param1, dealias = 1), init, param1 ) 

param2 = ( ϵ = 1/4, μ = Inf, N = 2^9, L = 2*π, T = 5, dt = 0.001 )

problem2 = Problem( WWn(param2, dealias = 1), init, param2 ) 

solve!([problem1, problem2])

plot(problem1; fourier = true)
plot!(problem2; fourier = true)
```

```@example fourier
plot([problem1, problem2]; fourier = true, label = :none)
```



