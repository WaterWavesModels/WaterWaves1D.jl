# # Workspace
#
#md # [`notebook`](@__NBVIEWER_ROOT_URL__notebooks/two_problems.ipynb)
#
#using ShallowWaterModels
include("../src/dependencies.jl")


#----
param = ( μ  = .1,
			ϵ  = .1,
        	N  = 2^11,
            L  = 20,
			c  = 1.05,
			T  = 80/1.05,
            dt = 0.001/1.05,
			q = 0.5
			)
mesh = Mesh(param)

function sol(x,c,ϵ,μ)
	(c^2-1)/(ϵ*c)*sech.(sqrt(3*(c^2-1)/μ)/2*x).^2
end

(η,v,flag) = SolitaryWaveWhithamBoussinesq(mesh, merge(param,(α  = 0,)), sol(mesh.x,c,ϵ,μ); iterative = true, q=1)
(η,v,flag) = SolitaryWaveWhithamBoussinesq(mesh, merge(param,(α  = 1/2,)), v; iterative = false, q=1)

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
