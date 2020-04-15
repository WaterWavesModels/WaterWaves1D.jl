"""
A proof of concept: numerically computes the solitary wave of the KdV equation with prescribed velocity
"""
#using ShallowWaterModels
include("../src/dependencies.jl")


#----
function test()
	param = ( μ  = 1,
			ϵ  = 1,
        	N  = 2^10,
            L  = 10,
			)
	mesh = Mesh(param)

	function sol(x,α)
		2*α*sech.(sqrt(3*2*α)/2*x).^2
	end

	u₀ = SolitaryWaveWhitham(mesh, merge(param,(c=2.5,)), sol(mesh.x,1.5)+.1*exp.(-(mesh.x).^2); iterative = true, KdV = true, α=1,tol = 1e-8)
	u₁ = SolitaryWaveWhitham(mesh, merge(param,(c=2.6,)), sol(mesh.x,1.6)+.1*exp.(-(mesh.x).^2); iterative = false, KdV = true, α=1,tol = 1e-8)
	u₂ = SolitaryWaveWhitham(mesh, merge(param,(c=2.7,)), sol(mesh.x,1.7)+.5*exp.(-(mesh.x).^2); iterative = false, KdV = true, α=1,tol = 1e-8, q=2/3)

	# barycentric Lagrange inerpolation (4.2) in https://people.maths.ox.ac.uk/trefethen/barycentric.pdf
	function P(X,dc,u₀,u₁,u₂)
		(u₀/(2*(X/dc+2))-u₁/(X/dc+1)+u₂/(2*X/dc))/(1/(2*(X/dc+2))-1/(X/dc+1)+1/(2*X/dc))
	end


	function solve!(u₀,u₁,u₂,c,dc)
			tu₂ = copy(u₂)
			u₂ .= SolitaryWaveWhitham(mesh, merge(param,(c=c,)), P.(dc,dc,u₀,u₁,u₂); iterative = false, α=1, KdV = true, verbose = true, max_iter = 5, tol = 1e-10)
			u₀ .= u₁
			u₁ .= tu₂
	end

	for cs = range(2.8; step = 0.1, stop = 5)
		print(string("c = ",cs,"\n"))
		solve!(u₀,u₁,u₂,cs,0.1)
	end

	ucomp = u₂
	uexact = sol(mesh.x,4)

	plt = plot(layout=(2,2))

    plot!(plt[1,1], mesh.x, [ucomp uexact];
     title=string("c=5"),
     label=["computed" "formula"])

    plot!(plt[2,1], fftshift(mesh.k),
     [log10.(abs.(fftshift(fft(ucomp)))) log10.(abs.(fftshift(fft(uexact))))];
     title="frequency",
	 label=["computed" "formula"])

	 plot!(plt[1,2], mesh.x, ucomp-uexact;
      title=string("c=5"),
      label="difference")

     plot!(plt[2,2], fftshift(mesh.k),
      log10.(abs.(fftshift(fft(ucomp-uexact))));
      title="frequency",
	  label="difference")
end
