# # Workspace
#
#md # [`notebook`](@__NBVIEWER_ROOT_URL__notebooks/two_problems.ipynb)
#
#using ShallowWaterModels
include("../src/dependencies.jl")


#----
param = ( μ  = .01,
			ϵ  = .01,
        	N  = 2^11,
            L  = 20,
			α  = 1/2,
			c  = 1.005,
			T  = 100,
            dt = 0.01,
			)
mesh = Mesh(param)

function sol(x,c,ϵ,μ)
	(c^2-1)/(ϵ*c)*sech.(sqrt(3*(c^2-1)/μ)/2*x).^2
end

(η,v,flag) = SolitaryWaveWhithamBoussinesq(mesh, param, sol(mesh.x,c,ϵ,μ); iterative = true)

init=Init(mesh,η,v)

model  = fdBoussinesq(param)
solver   = RK4(param,model)
problem = Problem(model, init, param, solver);

@time solve!( problem )

p = plot(layout=(2,1))
fig_problem!( p, problem )
display(p)

@time create_animation( problem )
