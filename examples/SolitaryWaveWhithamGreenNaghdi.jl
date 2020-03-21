# # Workspace
#
#md # [`notebook`](@__NBVIEWER_ROOT_URL__notebooks/two_problems.ipynb)
#
#using ShallowWaterModels
include("../src/dependencies.jl")


#----
param = ( μ  = 1,
			ϵ  = 1,
        	N  = 2^11,
            L  = 20,
			T  = 80/2,
            dt = 0.001/2
						)
mesh = Mesh(param)

function sol(x,c,ϵ,μ)
	(c^2-1)/ϵ*sech.(sqrt(3*(c^2-1)/(c^2)/μ)/2*x).^2
end
ηGN(c)= sol(mesh.x,c,param.ϵ,param.μ)
(η₀,v₀,flag) = SolitaryWaveWhithamGreenNaghdi(
				mesh, merge(param,(c=2,)), ηGN(2);
				tol =  1e-12, max_iter=20,
				iterative = true,
				verbose = true, )

(η₁,v₁,flag) = SolitaryWaveWhithamGreenNaghdi(
				mesh, merge(param,(c=3,)), ηGN(3);
				tol =  1e-12, max_iter=20,
				iterative = true,
				verbose = true, )

(η₂,v₂,flag) = SolitaryWaveWhithamGreenNaghdi(
				mesh, merge(param,(c=4,)), ηGN(4);
				tol =  1e-12, max_iter=20,
				iterative = true,
				verbose = true, )

# barycentric Lagrange inerpolation (4.2) in https://people.maths.ox.ac.uk/trefethen/barycentric.pdf
function P(X,dc,u₀,u₁,u₂)
	(u₀/(2*(X/dc+2))-u₁/(X/dc+1)+u₂/(2*X/dc))/(1/(2*(X/dc+2))-1/(X/dc+1)+1/(2*X/dc))
end

function Q(X,dc,u₀,u₁,u₂)
	(u₀./(2*(X/dc+2)) .-u₁./(X/dc+1) .+u₂./(2*X/dc))./(1/(2*(X/dc+2))-1/(X/dc+1)+1/(2*X/dc))
end

function solve2!(η₀,η₁,η₂,c,dc)
		tη₂ = copy(η₂)
		η₂ .= SolitaryWaveWhithamGreenNaghdi(mesh, merge(param,(c=c,)), P.(dc,dc,η₀,η₁,η₂);
				tol =  1e-10, max_iter=20,
				iterative = true, verbose = true)[1]
		η₀ .= η₁
		η₁ .= tη₂
end

function solve2!(η₀,η₁,η₂,c,old_dc,new_dc)
		tη₀ = copy(η₀)
		tη₁ = copy(η₁)
		tη₂ = copy(η₂)
		η₀ .= tu₂
		η₁ .= SolitaryWaveWhitham(mesh, merge(param,(c=c+1*new_dc,)), P.(1*new_dc,old_dc,tη₀,tη₁,tη₂);
			tol =  1e-10, max_iter=20,
			iterative = true, verbose = true)[1]
		η₂ .= SolitaryWaveWhitham(mesh, merge(param,(c=c+2*new_dc,)), P.(2*new_dc,old_dc,tη₀,tη₁,tη₂);
			tol =  1e-10, max_iter=20,
			iterative = true, verbose = true)[1]
end

for cs = range(5; step = 1, stop = 20)
	print(string("c = ",cs,"\n"))
	solve2!(η₀,η₁,η₂,cs,1)
end

solve2!(u₀,u₁,u₂,1.221,0.01,1e-4)


for cs = range(1.221+3e-4; step = 1e-4, stop = 1.229)
	print(string("c = ",cs,"\n"))
	solve2!(u₀,u₁,u₂,cs,1e-4)
end


plot(mesh.k,log10.(abs.(real.(fft(η₂)))))


init=Init(mesh,η,v)

model  = fdBoussinesq(param)
solver   = RK4(param,model)
problem = Problem(model, init, param, solver);

@time solve!( problem )

p = plot(layout=(2,1))
fig_problem!( p, problem )
display(p)

plot(mesh.x,[η-mapfro(model,problem.data.U[end])[1]])

@time create_animation( problem )
