export figure1,figure2,figure3,figure4,figure6
"""
Reproduces the figures in the work of C. Klein and V. Duchêne
"""
#using ShallowWaterModels
include("../src/dependencies.jl")


#---- Figure 1
"""
	figure1()

WGN solitary wave of velocity c = 1.1
"""
function figure1()
	param = ( μ  = 1,
			ϵ  = 1,
        	N  = 2^8,
            L  = 10*π,
						)
	mesh = Mesh(param)
	function sol(x,c,ϵ,μ)
		(c^2-1)/ϵ*sech.(sqrt(3*(c^2-1)/(c^2)/μ)/2*x).^2
	end
	ηGN(c)= sol(mesh.x,c,param.ϵ,param.μ)
	uGN(c)= c*ηGN(c)./(1 .+ param.ϵ*ηGN(c))

	(η,u) = SolitaryWaveWhithamGreenNaghdi(
				mesh, merge(param,(c=1.1,)), ηGN(1.1);
				method=2, α = 0,
				tol =  1e-16, max_iter=4,
				ktol =1e-11, gtol = 1e-16,
				iterative = true, q=1,
				verbose = false, SGN = false)

	plt = plot(layout=(1,2))
	plot!(plt[1,1], mesh.x, [u uGN(1.1)];
	  title="c=1.1",
	  label=["WGN" "SGN"])

	plot!(plt[1,2], fftshift(mesh.k),
	  [log10.(abs.(fftshift(fft(u))))
	  log10.(abs.(fftshift(fft(uGN(1.1)))))];
	  title="frequency",
	  label="WGN")
	display(plt)
	savefig("fig1.pdf");
end

#---- Figure 2
"""
	figure2(c)

WGN solitary wave of velocity `c` (optional, default is `c = 2`)
Uses GMRES-iterative method for the inversion of the Jacobian.
"""
function figure2(c=2)
	param = ( μ  = 1,
		ϵ  = 1,
      	N  = 2^10,
        L  = 10*π,)
	mesh = Mesh(param)
	function sol(x,c,ϵ,μ)
		(c^2-1)/ϵ*sech.(sqrt(3*(c^2-1)/(c^2)/μ)/2*x).^2
	end
	ηGN(c)= sol(mesh.x,c,param.ϵ,param.μ)
	uGN(c)= c*ηGN(c)./(1 .+ param.ϵ*ηGN(c))

	(η,u) = SolitaryWaveWhithamGreenNaghdi(
				mesh, merge(param,(c=c,)), ηGN(c);
				method=2, α = 0,
				tol =  1e-12, max_iter=10,
				iterative = true, q=1,
				verbose = false, SGN = false)
	plt = plot(layout=(1,2))
	plot!(plt[1,1], mesh.x, [u uGN(c)];
	  title=string("c=",c),
	  label=["WGN" "SGN"])
	plot!(plt[1,2], fftshift(mesh.k),
	  [log10.(abs.(fftshift(fft(u)))) log10.(abs.(fftshift(fft(uGN(c)))))];
	  title="frequency",
	  label=["WGN" "SGN"])
	display(plt)
	savefig("fig2.pdf");
end

#---- Figure 3
"""
	figure3(c)

WGN solitary wave of velocity `c` (optional, default is `c = 20`)
Uses non-iterative method for the inversion of the Jacobian.
"""
function figure3(c=20)
	param = ( μ  = 1,
		ϵ  = 1,
      	N  = 2^10,
        L  = 10*π,)
	mesh = Mesh(param)
	function sol(x,c,ϵ,μ)
		(c^2-1)/ϵ*sech.(sqrt(3*(c^2-1)/(c^2)/μ)/2*x).^2
	end
	ηGN(c)= sol(mesh.x,c,param.ϵ,param.μ)
	uGN(c)= c*ηGN(c)./(1 .+ param.ϵ*ηGN(c))

	(η,u) = SolitaryWaveWhithamGreenNaghdi(
				mesh, merge(param,(c=c,)), ηGN(c);
				method=2, α = 1, #ici α = 1 évite des oscillations important si c = 3 ou c = 20
				tol =  1e-14, max_iter=15,
				iterative = false, q=1,
				verbose = true, SGN = false)

	plt = plot(layout=(1,2))
	plot!(plt[1,1], mesh.x, [u uGN(c)];
	  title=string("c=",c),
	  label=["WGN" "SGN"])
	plot!(plt[1,2], fftshift(mesh.k),
	  [log10.(abs.(fftshift(fft(u)))) log10.(abs.(fftshift(fft(uGN(c)))))];
	  title="frequency",
	  label=["WGN" "SGN"])
    display(plt)
	savefig("fig3.pdf");
end
#---- Figure 4
"""
	figure4(c)

WGN solitary wave of velocity `c` (optional, default is `c = 100`)
Uses non-iterative method for the inversion of the Jacobian.
Uses the rescaled variable equation.
"""
function figure4(c=100)
	param = ( μ  = 1,
		ϵ  = 1,
    	N  = 2^10,
      	L  = 10*π,)
	mesh = Mesh(param)
	function sol(x,c,ϵ,μ)
		(c^2-1)/ϵ*sech.(sqrt(3*(c^2-1)/(c^2)/μ)/2*x).^2
	end
	ηGN(c)= sol(mesh.x,c,param.ϵ,param.μ)
	uGN(c)= c*ηGN(c)./(1 .+ param.ϵ*ηGN(c))

	(η,u) = SolitaryWaveWhithamGreenNaghdi(
				mesh, merge(param,(c=c,)), ηGN(c);
				method=3, α = 1,
				tol =  1e-10, max_iter=10,
				iterative = false, q=1,
				verbose = true, SGN = false)

	plt = plot(layout=(1,2))
	plot!(plt[1,1], mesh.x, [u/c uGN(c)/c];
	  title=string("c=",c),
	  label=["WGN" "SGN"])
	plot!(plt[1,2], fftshift(mesh.k),
	  [log10.(abs.(fftshift(fft(u/c))))
	  log10.(abs.(fftshift(fft(uGN(c)/c))))];
	  title="frequency",
	  label="WGN")

	savefig("fig4.pdf");
end

#------ Figure 6
"""
	figure6()

Jacobian associated with WGN solitary wave of velocity 20.
"""
function figure6()
	c=20
	ϵ,μ,α=1,1,0
	L,N=10*π,2^10
	mesh = Mesh(L,N)
	k,x=mesh.k,mesh.x
	F₁ 	= tanh.(sqrt(μ)*abs.(k))./(sqrt(μ)*abs.(k))
	F₁[1] 	= 1
	F₀       = 1im * sqrt.(3*(1 ./F₁ .- 1)).*sign.(k)
	x₀ = mesh.x[1]
	FFT = exp.(-1im*k*(x.-x₀)');
	IFFT = exp.(1im*k*(x.-x₀)')/length(x);
	M₀ = IFFT * Diagonal( F₀ )* FFT
	M(v) = Diagonal( v )
	function sol(x,c,ϵ,μ)
		(c^2-1)/ϵ*sech.(sqrt(3*(c^2-1)/(c^2)/μ)/2*x).^2
	end
	guess = sol(mesh.x,c,ϵ,μ)
	u = c*guess./(1 .+ ϵ*guess)
	dxu = real.(ifft(1im*k.*fft(u)))
	dxu ./= norm(dxu,2)
	hu = c ./(c .- ϵ*u)
	Fu = hu.* real.(ifft(F₀.*fft(u)))
	F2u = real.(ifft(F₀.*fft(hu.^2 .* Fu ) ))
	Du = c ./hu .+ 2*ϵ/3 * F2u ./ hu .+ ϵ^2/c * hu .* Fu.^2 .- (hu.^2)/c
	Jac = (-1/3 *M(1 ./ hu.^2 )* M₀ * M(hu.^3)* (c*M₀ .+ 3*ϵ * M( Fu ))
					.+ ϵ * M( hu .* Fu ) * M₀
					.+ M( Du ) .+ α*M( hu.^2 )*dxu*dxu' *M(1 ./ hu.^2 ) )
	Jacstar = -1/3 *M(1 ./ hu.^2 )* M₀ * M(hu.^3)* c*M₀

	plt = plot(layout=(1,2))
	surface!(plt[1,1],fftshift(k),fftshift(k)[N:-1:1],log10.(abs.(FFT*Jac*IFFT)))
	surface!(plt[1,2],fftshift(k),fftshift(k)[N:-1:1],log10.(abs.(FFT*Jacstar*IFFT)))

	savefig("fig6.pdf");
end
