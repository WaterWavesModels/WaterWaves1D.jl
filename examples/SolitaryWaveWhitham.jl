# # Workspace
#
#md # [`notebook`](@__NBVIEWER_ROOT_URL__notebooks/two_problems.ipynb)
#
#using ShallowWaterModels
include("../src/dependencies.jl")


#----
param = ( μ  = 1,
			ϵ  = 1,
        	N  = 2^10,
            L  = 10,
			)
mesh = Mesh(param)

function sol(x,α)
	2*α*sech.(sqrt(3*2*α)/2*x).^2
end

u₀ = SolitaryWaveWhitham(mesh, merge(param,(c=1.1,)), sol(mesh.x,0.1); iterative = true)
u₁ = SolitaryWaveWhitham(mesh, merge(param,(c=1.101,)), sol(mesh.x,0.1); iterative = true)


for cs = range(1.101; step = 0.001, stop = 1.2)
	print(string("c = ",cs,"\n"))
	u₀ .= u₁
	u₁ .= SolitaryWaveWhitham(mesh, merge(param,(c=cs,)), u₁+u₁-u₀; iterative = true, verbose = false)
end

for cs = range(1.2001; step = 0.0001, stop = 1.225)
	print(string("c = ",cs,"\n"))
	u₀ .= u₁
	u₁ .= SolitaryWaveWhitham(mesh, merge(param,(c=cs,)), u₁+1/10*(u₁-u₀); iterative = true, verbose = false)
end

for cs = range(1.22501; step = 0.00001, stop = 1.229)
	print(string("c = ",cs,"\n"))
	u₀ .= u₁
	u₁ .= SolitaryWaveWhitham(mesh, merge(param,(c=cs,)), u₁+1/100*(u₁-u₀); iterative = true, verbose = false)
end

plot(mesh.x,u₁)
