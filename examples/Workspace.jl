# # Workspace
#
#md # [`notebook`](@__NBVIEWER_ROOT_URL__notebooks/two_problems.ipynb)
#
using ShallowWaterModels,Plots,FFTW
#include("../src/dependencies.jl")

#----

param = ( μ  = 1, ϵ  = 1, N  = 2^13,
		h₀=1.0962,h₁=1.1,h₂=1.2)
(η,u,v,mesh,par)=CnoidalWaveWhithamGreenNaghdi(param;P=160,SGN=true)
x=mesh.x
using Statistics
h=1 .+η
meanη=mean(h)
m=-sqrt(param.h₀*param.h₁*param.h₂)

M=-m*mean(1 ./ h)/sqrt(mean(h)) #Froude number, close to 1

L=80*par.λ

interp(x)=sech.(x/3).^2
η₀(x)=(0.13173.+(1.2-1.13173)*interp(x))
η1=η.*(x.>=-L).*(x.<=L)+η₀.(x.+L).*(x.<-L)+η₀.(x.-L).*(x.>L)
v₀(x)=(meanv .+(v[Int(end/2)]-meanv)*interp(x))
v1=v.*(x.>=-L).*(x.<=L)+v₀(x.+L).*(x.<-L)+v₀(x.-L).*(x.>L)
v1 .-=par.c


interp(x)=sech.(x/2).^2
η₀(x)=(0.13173.+(1.2-1.13173)*interp(x))
η2=η.*(x.>=-L).*(x.<=L)+η₀.(x.+L).*(x.<-L)+η₀.(x.-L).*(x.>L)
v₀(x)=(meanv .+(v[Int(end/2)]-meanv)*interp(x))
v2=v.*(x.>=-L).*(x.<=L)+v₀(x.+L).*(x.<-L)+v₀(x.-L).*(x.>L)
v2 .-=par.c


interp(x)=sech.(10*x).^2
η₀(x)=(0.13173.+(1.2-1.13173)*interp(x))
η3=η.*(x.>=-L).*(x.<=L)+η₀.(x.+L).*(x.<-L)+η₀.(x.-L).*(x.>L)
v₀(x)=(meanv .+(v[Int(end/2)]-meanv)*interp(x))
v3=v.*(x.>=-L).*(x.<=L)+v₀(x.+L).*(x.<-L)+v₀(x.-L).*(x.>L)
v3 .-=par.c


plot(mesh.x,[η1 η2 η3])
f1=fftshift(log10.(abs.(fft(η1))))
f2=fftshift(log10.(abs.(fft(η2))))
f3=fftshift(log10.(abs.(fft(η3))))

plot(fftshift(mesh.k),[f1 f2 f3])

init1=Init(mesh,η1,v1)
init2=Init(mesh,η2,v2)
init3=Init(mesh,η3,v3)
param2=merge(param,(L=-mesh.xmin,T  = 100, dt = 0.01,))
model=WhithamGreenNaghdi(param2,SGN=true)
problem1 = Problem(model, init1, param2)
problem2 = Problem(model, init2, param2)
problem3 = Problem(model, init3, param2)
solve!( problem1 )
solve!( problem2 )
solve!( problem3 )

p = plot(layout=(2,1))
fig_problem!(p,problem1,50)
fig_problem!(p,problem2,50)
fig_problem!(p,problem3,50)
savefig("Gavrilyuk50.pdf")
xlims!(-L-100,-L+10)
savefig("Gavrilyuk50zoom.pdf")


p = plot(layout=(2,1))
fig_problem!(p,problem1,100)
fig_problem!(p,problem2,100)
fig_problem!(p,problem3,100)
savefig("Gavrilyuk100.pdf")
xlims!(-L-100,-L+10)
savefig("Gavrilyuk100zoom.pdf")

save(problem1,"problem1")
save(problem2,"problem2")
save(problem3,"problem3")


#----
param = ( μ  = 0.5,
			ϵ  = 0.5,
        	N  = 2^9,
            L  = 10,
            T  = 5,
            dt = 0.01,
			α = 1/2,)

init     = BellCurve(merge(param,(θ = 2,)))
model = WhithamGreenNaghdi(param;iterate=true)
problem1 = Problem(model, init, param)

@time solve!( problem1 )

model2 = fdBoussinesq(param)
problem2 = Problem(model2, init, param)
@time solve!( problem2 )

model3 = WaterWaves(param)
problem3 = Problem(model3, init, param)
@time solve!( problem3 )

p = plot(layout=(2,1))
fig_problem!( p, problem1, 10)
fig_problem!( p, problem2, 10)
fig_problem!( p, problem3, 10)

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
