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

c=1.101

dc = 0.01

function solve!(u₀,u₁,c,dc)
	while dc > 0.00001
		new_dc = dc
		@info string("c = ",c,"dc = ",dc,"\n")
		u₀[1] .= u₁[1]
		flag = 1
		while flag==1
			u₁ = SolitaryWaveWhitham(mesh, merge(param,(c=c+new_dc,)), u₁[1]+new_dc/dc*(u₁[1]-u₀[1]); iterative = true, verbose = false)
			if u₁[2]==1
				new_dc = new_dc/2
			else
				flag=0
			end
		end
		c += new_dc
		dc = new_dc
	end
	return (u₀,u₁,c,dc)
end

#(u₀,u₁,c,dc)=solve!(u₀,u₁,c,dc)

function solve!(u₀,u₁,c,old_dc,new_dc)
		u₀[1] .= u₁[1]
		u₁[1] .= SolitaryWaveWhitham(mesh, merge(param,(c=c,)), u₁[1]+new_dc/old_dc*(u₁[1]-u₀[1]); iterative = true, verbose = false)[1]
end

for cs = range(1.101; step = 0.001, stop = 1.2)
	print(string("c = ",cs,"\n"))
	solve!(u₀,u₁,cs,0.001,0.001)
end

for cs = range(1.2; step = 0.0005, stop = 1.225)
	print(string("c = ",cs,"\n"))
	solve!(u₀,u₁,cs,0.0001,0.0001)
end

for cs = range(1.225; step = 0.00001, stop = 1.229)
	print(string("c = ",cs,"\n"))
	solve!(u₀,u₁,cs,0.0001,0.0001)
end

for cs = range(1.229; step = 0.000001, stop = 1.22904)
	print(string("c = ",cs,"\n"))
	solve!(u₀,u₁,cs,0.00001,0.00001)
end

plot(mesh.x,u₁[1])
plot(mesh.k,log10.(abs.(real.(fft(u₁[1])))))
