# #
# A proof of concept:
# 1) PlotSolitaryWaveWhithamBoussinesq() numerically computes the solitary wave of a Whitham-Boussinesq equation with prescribed velocity
# 2) IntegrateSolitaryWaveWhithamBoussinesq() integrates the equation in time
# #
export PlotSolitaryWaveWhithamBoussinesq,IntegrateSolitaryWaveWhithamBoussinesq
using WaterWaves1D,FFTW,Plots;
@info("Construct solitary wave solutions of a Whitham-Boussinesq equation (use `PlotSolitaryWaveWhithamBoussinesq`).
Integrate in time the equations starting from such initial data (use  `IntegrateSolitaryWaveWhithamBoussinesq`.
This allows to assess the validity and precision of the numerical scheme.")

#----
"""
	PlotSolitaryWaveWhithamBoussinesq(;kwargs...)

Construct and plot a solitary wave of the Whitham-Boussinesq equation.

All keyword arguments are optional.
- `c`, the velocity of the wave (default `c=1.05`),
- `α`, parameter of the Whitham-Boussinesq model (default `α=1`),
- `L`, the half-length of the mesh (default `L=20`),
- `N`, the number of collocation points (default `N=2^9`),
- `μ`, the shallowness parameter (default `μ=0.1`),
- `ϵ`, the nonlinearity parameter (default `ϵ=0.1`).

Return the difference between the computed solution and the solution of the Serre-Green-Naghdi equation (in ℓ^∞ norm at collocation points).
"""
function PlotSolitaryWaveWhithamBoussinesq(;c=1.05,α=1,L=20,N=2^9,μ=0.1,ϵ=0.1)
	param = ( μ  = μ, ϵ = ϵ,
			N = N, L = L,
			c = c)

	(η,v,mesh) = SolitaryWaveWhithamBoussinesq( param;
										β = 1, α  = α, iterative = false, max_iter=30)

	ηSGN = (c^2-1)/ϵ*sech.(sqrt(3*(c^2-1)/(c^2)/μ)/2*(mesh.x)).^2
	uSGN = c*η./(1 .+ ϵ*η)
	k = mesh.k
	F₀ = sqrt(param.μ)*1im * k
	DxF(v) = real.(ifft(F₀ .* fft(v)))
	h = 1 .+ param.ϵ*ηSGN
	vSGN = uSGN - 1/3 ./h .* (DxF(h.^3 .*DxF(uSGN)))


	plt = plot(layout=(1,2))
	plot!(plt[1,1], mesh.x, [η v];
    	title=string("initial data"),
    	label=["η" "v"])
	plot!(plt[1,1], mesh.x, [ηSGN vSGN];
    	title=string("initial data"),
      	label=["ηSGN" "vSGN"])


    plot!(plt[1,2], fftshift(mesh.k),
    	[abs.(fftshift(fft(η))).+eps() abs.(fftshift(fft(v))).+eps()];
     	yscale=:log10,
	 	title="frequency",
	 	label=["|η̂|" "|v̂|"])
	plot!(plt[1,2], fftshift(mesh.k),
      	[abs.(fftshift(fft(ηSGN))).+eps() abs.(fftshift(fft(vSGN))).+eps()];
	  	yscale=:log10,
	  	title="frequency",
 	 	label=["|η̂SGN|" "|v̂SGN|"])

	display(plt)
	return 	maximum(abs.(η-ηSGN))

end

"""
	IntegrateSolitaryWaveWhithamBoussinesq(;kwargs...)

Integrate in time a Whitham-Boussinesq equation with a solitary wave initial data

All keyword arguments are optional.
- `c`, the velocity of the wave (default `c=1.05`),
- `α`, parameter of the Whitham-Boussinesq model (default `α=1`),
- `L`, the half-length of the mesh (default `L=20`),
- `T`, the final computation time (default `T=40/1.05`),
- `dt`, the timestep (default `T=0.001/1.05`),
- `ns`, storing data every `ns` computation steps (default `ns=200`),
- `anim`, builds an animation if `true` (default),
- `N`, the number of collocation points (default `N=2^9`),
- `μ`, the shallowness parameter (default `μ=0.1`),
- `ϵ`, the nonlinearity parameter (default `ϵ=0.1`).

Return the difference between the initial data and final data (in ℓ^∞ norm at collocation points).
"""
function IntegrateSolitaryWaveWhithamBoussinesq(; 
	α  = 1,
	μ  = 1,
	ϵ  = 1,
	c = 1.05,
	N  = 2^9,
	L  = 20,
	T  = 40/1.05,
	dt = 0.001/1.05,
	ns = 200,
	anim = true)
	param = ( μ  = μ,
			ϵ  = ϵ,
			c = c,
        	N  = N,
            L  = L,
			T  = T,
            dt = dt,
			ns = ns
			)

	(η,v,mesh) = SolitaryWaveWhithamBoussinesq( param;
										β = 1, α  = α, iterative = false, max_iter=30)

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

	model  = WhithamBoussinesq(param;α=α)
	problem = Problem(model, init, param);

	solve!( problem )
	ηf, = solution(problem)
	p = plot(layout=(1,2))
	plot!(p[1,1],mesh.x,[η ηf];
		title="surface deformation",
		label=["initial" "final"])
	plot!(p[1,2],mesh.x,η-ηf;
		title="difference",
		label="")
	display(p)
	if anim ==true
		@gif for t in problem.times.ts
			plot(problem, T = t)
		end
	end
	
	return 	maximum(abs.(ηf-η))


end
nothing
