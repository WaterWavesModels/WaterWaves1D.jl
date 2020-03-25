"""
Reproduces the figures in the work of C. Klein and V. Duchêne
"""
#using ShallowWaterModels
include("../src/dependencies.jl")

function sol(x,c,ϵ,μ)
	(c^2-1)/ϵ*sech.(sqrt(3*(c^2-1)/(c^2)/μ)/2*x).^2
end
ηGN(c)= sol(mesh.x,c,param.ϵ,param.μ)
uGN(c)= c*ηGN(c)./(1 .+ param.ϵ*ηGN(c))


#---- Figure 1
function figure1()
	param = ( μ  = 1,
			ϵ  = 1,
        	N  = 2^8,
            L  = 10*π,
						)
	mesh = Mesh(param)
	ηGN(c)= sol(mesh.x,c,param.ϵ,param.μ)
	uGN(c)= c*ηGN(c)./(1 .+ param.ϵ*ηGN(c))

	(η,u) = SolitaryWaveWhithamGreenNaghdi(
				mesh, merge(param,(c=1.1,)), ηGN(1.1);
				method=2, α = 0,
				tol =  1e-16, max_iter=4,
				ktol =1e-11, gtol = 1e-14,
				iterative = true, q=1,
				verbose = false, GN = false)

	plt = plot(layout=(1,2))
	plot!(plt[1,1], mesh.x, [u uGN(1.1)];
	  title="c=1.1",
	  label=["WGN" "SGN"])

	plot!(plt[1,2], fftshift(mesh.k),
	  log10.(abs.(fftshift(fft(u))));
	  title="frequency",
	  label="WGN")
end

#---- Figure 2
function figure2()
	param = ( μ  = 1,
		ϵ  = 1,
      	N  = 2^10,
        L  = 10*π,)
	mesh = Mesh(param)
	ηGN(c)= sol(mesh.x,c,param.ϵ,param.μ)
	uGN(c)= c*ηGN(c)./(1 .+ param.ϵ*ηGN(c))

	(η,u) = SolitaryWaveWhithamGreenNaghdi(
				mesh, merge(param,(c=2,)), ηGN(2);
				method=2, α = 0,
				tol =  1e-10, max_iter=10,
				iterative = false, q=1,
				verbose = true, GN = false)
	plt = plot(layout=(1,2))
	plot!(plt[1,1], mesh.x, [u uGN(2)];
	  title="c=2",
	  label=["WGN" "SGN"])
	plot!(plt[1,2], fftshift(mesh.k),
	  log10.(abs.(fftshift(fft(u))));
	  title="frequency",
	  label="WGN")
end

#---- Figure 3
function figure3()
	param = ( μ  = 1,
		ϵ  = 1,
      	N  = 2^10,
        L  = 10*π,)
	mesh = Mesh(param)
	ηGN(c)= sol(mesh.x,c,param.ϵ,param.μ)
	uGN(c)= c*ηGN(c)./(1 .+ param.ϵ*ηGN(c))

	(η,u) = SolitaryWaveWhithamGreenNaghdi(
				mesh, merge(param,(c=20,)), ηGN(20);
				method=2, α = 1,
				tol =  1e-10, max_iter=15,
				iterative = false, q=1,
				verbose = true, GN = false)

	plt = plot(layout=(1,2))
	plot!(plt[1,1], mesh.x, [u uGN(20)];
	  title="c=20",
	  label=["WGN" "SGN"])
	plot!(plt[1,2], fftshift(mesh.k),
	  log10.(abs.(fftshift(fft(u))));
	  title="frequency",
	  label="WGN")
end
#---- Figure 4
function figure4()
	param = ( μ  = 1,
		ϵ  = 1,
    	N  = 2^10,
      	L  = 10*π,)
	mesh = Mesh(param)
	ηGN(c)= sol(mesh.x,c,param.ϵ,param.μ)
	uGN(c)= c*ηGN(c)./(1 .+ param.ϵ*ηGN(c))

	(η,u) = SolitaryWaveWhithamGreenNaghdi(
				mesh, merge(param,(c=100,)), ηGN(100);
				method=3, α = 0,
				tol =  1e-10, max_iter=40,
				iterative = false, q=1,
				verbose = true, GN = false)

	plt = plot(layout=(1,2))
	plot!(plt[1,1], mesh.x, [u/100 uGN(100)/100];
	  title="c=100",
	  label=["WGN" "SGN"])
	plot!(plt[1,2], fftshift(mesh.k),
	  log10.(abs.(fftshift(fft(u/100))));
	  title="frequency",
	  label="WGN")
  end

#------ Figure 5
function figure5()
	c=20
	ϵ,μ,α=1,1,1
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
	guess = sol(mesh.x,c,param.ϵ,param.μ)
	u = c*guess./(1 .+ ϵ*guess)
	dxu = real.(ifft(1im*k.*fft(u)))
	dxu ./= norm(dxu,2)
	hu = c ./(c .- ϵ*u)
	Fu = hu.* real.(ifft(F₀.*fft(u)))
	F2u = real.(ifft(F₀.*fft(hu.^2 .* Fu ) ))
	Du = c ./hu .+ 2*ϵ/3 * F2u ./ hu .+ ϵ^2/c * hu .* Fu.^2 .- (hu.^2)/c
	Jac = (-1/3 *M(1 ./ hu.^2 )* M₀ * M(hu.^3)* (c*M₀ .+ 3*ϵ * M( Fu ))
					.+ ϵ * M( hu .* Fu ) * M₀
					.+ M( Du ) .+ α*dxu*dxu' *M(1 ./ hu.^2 ) )
	Jacstar = -1/3 *M(1 ./ hu.^2 )* M₀ * M(hu.^3)* c*M₀

	plt = plot(layout=(1,2))
	surface!(plt[1,1],fftshift(k),fftshift(k)[N:-1:1],log10.(abs.(FFT*Jac*IFFT)))
	surface!(plt[1,2],fftshift(k),fftshift(k)[N:-1:1],log10.(abs.(FFT*Jacstar*IFFT)))
	return Jac
end

#------ Other experiments
# barycentric Lagrangian interpolation (4.2) in https://people.maths.ox.ac.uk/trefethen/barycentric.pdf
function P(X,dc,u₀,u₁,u₂)
	(u₀/(2*(X/dc+2))-u₁/(X/dc+1)+u₂/(2*X/dc))/(1/(2*(X/dc+2))-1/(X/dc+1)+1/(2*X/dc))
end

# function solve!(η₀,η₁,η₂,c,old_dc,new_dc)
# 		tη₀ = copy(η₀)
# 		tη₁ = copy(η₁)
# 		tη₂ = copy(η₂)
# 		η₀ .= tη₂
# 		η₁ .= SolitaryWaveWhithamGreenNaghdi(mesh, merge(param,(c=c+1*new_dc,)), P.(1*new_dc,old_dc,tη₀,tη₁,tη₂);
# 			tol =  1e-12, max_iter=20,
# 			iterative = false, verbose = true)[1]
# 		η₂ .= SolitaryWaveWhithamGreenNaghdi(mesh, merge(param,(c=c+2*new_dc,)), P.(2*new_dc,old_dc,tη₀,tη₁,tη₂);
# 			tol =  1e-12, max_iter=20,
# 			iterative = false, verbose = true)[1]
# end

function solve(c,method)
	SolitaryWaveWhithamGreenNaghdi(
				mesh, merge(param,(c=c,)), ηGN(c),
				method = method ;
				tol =  1e-10, max_iter=10,
				iterative = false,
				verbose = true, GN = false)[1]
end


function solve!(η₀,η₁,η₂,c,dc,method)
		tη₂ = copy(η₂)
		η₂ .= SolitaryWaveWhithamGreenNaghdi(mesh, merge(param,(c=c,)), P.(dc,dc,η₀,η₁,η₂), method=method;
				tol =  1e-10, max_iter=5, q=1, α=1,
				iterative = false, verbose = false, GN = false)[1]
		η₀ .= η₁
		η₁ .= tη₂
end

#------- Calculs
η₀ = solve(17,3)
η₁ = solve(18,3)
η₂ = solve(19,3)
for cs = range(20; step = 1, stop = 100)
	print(string("c = ",cs,"\n"))
	solve!(η₀,η₁,η₂,cs,1,1)
end

plot(mesh.x,[η η₂	])

plt = plot(layout=(1,2))
plot!(plt[1,1], mesh.x, [sη₂-η₂	])

plot!(plt[1,2], fftshift(mesh.k),
	  log10.(abs.(fftshift(fft(η₂))))
	  )

sη₁=copy(η₁);
sη₂=copy(η₂);
sη₃=copy(η₃);

η₀ = real.(ifft(interpolate(fft(sη₁),mesh.k)))




plot(mesh.x,η₀)

plot(mesh.k,log10.(abs.(real.(fft(η₁)))))

function interpolate(η,k)
	z=Complex.(zeros(length(k)))
	z[1:Int(length(η)/2)].=η[1:Int(length(η)/2)]
	z[length(z):-1:length(z)-Int(length(η)/2)+1].=η[length(η):-1:length(η)-Int(length(η)/2)+1]
	return length(k)/length(η).*z
end




#---------------- Evolution in time
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
