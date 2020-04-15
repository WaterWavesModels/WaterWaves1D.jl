"""
A proof of concept:
1) numerically computes the solitary wave of a Whitham-Boussinesq equation with prescribed velocity
2) integrates the equation with time
"""
#using ShallowWaterModels
include("../src/dependencies.jl")


#----
function test()
	param = ( μ  = .1,
			ϵ  = .1,
        	N  = 2^9,
            L  = 20,
			c  = 1.05,
			T  = 40/1.05,
            dt = 0.001/1.05,
			α	=1/2,
			)
	mesh = Mesh(param)

	function solKdV(x,c,ϵ,μ)
		(c^2-1)/(ϵ*c)*sech.(sqrt(3*(c^2-1)/μ)/2*x).^2
	end

	(η,v) = SolitaryWaveWhithamBoussinesq(mesh, param, solKdV(mesh.x,param.c,param.ϵ,param.μ)*param.c;
										α = 1, model  = 1/2, iterative = false, max_iter=30)

	init=Init(mesh,η,v)

	plt = plot(layout=(1,2))
	plot!(plt[1,1], mesh.x, [η v];
     title=string("initial data"),
     label=["η" "u"])

    plot!(plt[1,2], fftshift(mesh.k),
     [log10.(abs.(fftshift(fft(η)))) log10.(abs.(fftshift(fft(v))))];
     title="frequency",
	 label=["η" "u"])
	 display(plt)


	model  = fdBoussinesq(param)
	solver   = RK4(param,model)
	problem = Problem(model, init, param, solver);

	@time solve!( problem )
	ηf = mapfro(model,problem.data.U[end])[1]
	p = plot(layout=(1,2))
	plot!(p[1,1],mesh.x,[η ηf];
	title="surface deformation",
	label=["initial" "final"])
	plot!(p[1,2],mesh.x,η-ηf;
	title="difference",
	label="")
	display(p)


	#@time create_animation( problem )
end
