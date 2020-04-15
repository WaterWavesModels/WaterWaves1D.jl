"""
Minimal example for WhithamGreenNaghdi
"""
#using ShallowWaterModels
include("../src/dependencies.jl")

"""
	solveandplot()

Solves an initial value problem with WhithamGreenNaghdi model.
See parameters below.

"""
function solveandplot()
	# Set up
	param = ( μ  = 1, ϵ  = 1,
				N  = 2^10,   # number of collocation points
	            L  = 10*π,   # mesh of size 2L
	            T  = 1,      # final time
	            dt = 1/10^4, # timestep
				ns=1000      # stores data every ns time steps
				)

	mesh=Mesh(param); # construct mesh of collocation points, Fourier modes, etc.
	η= exp.(-(mesh.x).^2);v= 0*η; # Initial data

	init     = Init(mesh,η,v)
	model = WhithamGreenNaghdi(param;
				iterate=false) # iterate = true for GMRES iterate solver
	problem = Problem(model, init, param)

	# Compute
	@time solve!( problem )

	# Plot
	fftηfin=last(problem.data.U)[:,1]
	plt = plot(layout=(1,2))
	plot!(plt[1,1],fftshift(mesh.k),fftshift(log10.(abs.(fftηfin)));
			title = "frequencies")
	plot!(plt[1,2],mesh.x,real.(ifft(fftηfin));
			title = "physical space")
	display(plt)
end

solveandplot()
