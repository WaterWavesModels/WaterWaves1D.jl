# # Workspace
#
#md # [`notebook`](@__NBVIEWER_ROOT_URL__notebooks/two_problems.ipynb)
#
#using ShallowWaterModels
include("../src/dependencies.jl")

#----
param = ( μ  = .2,
			ϵ  = 0.2,
        	N  = 2^8,
            L  = 10,
            T  = 5,
            dt = 0.01,
			α = 1/2,
			SGN = false,)

init     = BellCurve(merge(param,(θ = 2,)))
model = WhithamGreenNaghdi(param)

problem = Problem(model, init, param)

@time solve!( problem )

model2 = fdBoussinesq(param)
problem2 = Problem(model2, init, param)
@time solve!( problem2 )

model3 = WaterWaves(param)
problem3 = Problem(model3, init, param)
@time solve!( problem3 )

p = plot(layout=(2,1))
fig_problem!( p, problem, 10)
display(p)


#----
param = ( μ  = 1,
			ϵ  = 1/1000,
        	N  = 2^9,
            L  = 10,
            T  = 10,
            dt = 0.001)
# k=Mesh(param).k
# x=Mesh(param).x
# ϵ  = 1/10
#
# 			f= k-> exp.(-k.*k)+ϵ*(abs.(ϵ^2*k).<1 ).*(rand(Float64,size(k)).-1/2)
# 			ζ=real(ifftshift(ifft(f(k))))
# 			plot(x,ζ)
# 			plot(k,f(k))

init     = BellCurve(merge(param,(θ = 2,)))
#init = Init((η = x-> exp.(-x.*x) .+ exp.(-100*x.*x)/10 , v = x->0 .*x))
#η = x -> 2 .^(-abs.((x).^4))
#v = x -> -10*x.*η(x)
#init = Init(η,v)


models = []
push!(models,PseudoSpectral(merge(param,(order=1,))))
push!(models,PseudoSpectral(merge(param,(order=2,))))
push!(models,PseudoSpectral(merge(param,(order=3,))))

#push!(models,Boussinesq(merge(param,(a=-1/3,b=1/3))))
push!(models,WaterWaves(param))

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
   	fig_problem!( p, problem, 10)
end
plot!(xlims=[(0,1) (-300,300 )])
plot!(ylims=[(-0*.005,0.005) (-15,5 )])
display(p)

savefig("final.pdf"); nothing # hide
#create_animations(problems)
#p=problems[1]
#plot(p.mesh.x,mapfro(p.model,p.data.U[100]))
