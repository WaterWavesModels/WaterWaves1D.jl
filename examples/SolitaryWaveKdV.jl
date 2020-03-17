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

u₀ = SolitaryWaveWhitham(mesh, merge(param,(c=2.5,)), sol(mesh.x,1.5)+1*exp.(-(mesh.x.-1).^2); iterative = true, KdV = true, tol = 1e-8)
u₁ = SolitaryWaveWhitham(mesh, merge(param,(c=2.6,)), sol(mesh.x,1.6); iterative = true, KdV = true, tol = 1e-8)

function solve!(u₀,u₁,c,old_dc,new_dc)
		u₀[1] .= u₁[1]
		u₁[1] .= SolitaryWaveWhitham(mesh, merge(param,(c=c,)), u₁[1]+new_dc/old_dc*(u₁[1]-u₀[1]); iterative = true, verbose = true, KdV = true, tol = 1e-10)[1]
end

for cs = range(2.6; step = 0.01, stop = 3)
	print(string("c = ",cs,"\n"))
	solve!(u₀,u₁,cs,0.01,0.01)
end

plot(mesh.x,[u₁[1] sol(mesh.x,2)])
plot(mesh.x,[u₁[1]-sol(mesh.x,2)])

plot(mesh.k,log10.(abs.(real.(fft(u₁[1])))))
plot(mesh.k,log10.(abs.(real.(fft(sol(mesh.x,1))))))
