# # Saved Experiments
#
#md # [`notebook`](@__NBVIEWER_ROOT_URL__notebooks/two_problems.ipynb)
#
using WaterWaves1D,Plots,FFTW,Statistics;gr()
#include("../src/dependencies.jl")

#---- Cold-Hot states transition for Serre-Green-Naghdi, discussed in the paper of Sergey Gavrilyuk, Boniface Nkonga, Keh-Ming Shyue and Lev Truskinovsky
function ColdHot(;scenario)
	if scenario !=1 && scenario !=2
		error("the scenario must be 1 or 2")
	end
	param = ( μ  = 1, ϵ  = 1, N  = 2^13,
			h₀=1.0962,h₁=1.1,h₂=1.2)
	@time (η,u,v,mesh,par)=CnoidalWaveWhithamGreenNaghdi(param;P=500,SGN=true)
	x=mesh.x
	h=1 .+ param.ϵ*η
	m=-sqrt(param.h₀*param.h₁*param.h₂)
	if scenario == 1
		meanη=mean(h)-1
		meanu=m*mean(1 ./h)
	else
		meanη=1.09808-1
		meanu=m/1.09808
	end

	M=-m*mean(1 ./ h)/sqrt(mean(h)) #Froude number, close to 1
	L=250*par.λ

	u.-=mean(u)
	u.+=m*mean(1 ./h)

	interp(x,δ)=sech.(x/δ).^2
	η₀(x,δ)=(meanη.+(0.2-meanη)*interp(x,δ))
	u₀(x,δ)=(meanu .+(u[Int(end/2)]-meanu)*interp(x,δ))

	Dx(v)=real.(ifft(1im*mesh.k.*fft(v)))

	η1=η.*(x.>=-L).*(x.<=L)+η₀(x.+L,0.1).*(x.<-L)+η₀(x.-L,0.1).*(x.>L)
	u1=u.*(x.>=-L).*(x.<=L)+u₀(x.+L,0.1).*(x.<-L)+u₀(x.-L,0.1).*(x.>L)
	v1=u1-param.μ/3*Dx(h.^3 .*Dx(u1))./h

	η2=η.*(x.>=-L).*(x.<=L)+η₀(x.+L,2).*(x.<-L)+η₀(x.-L,2).*(x.>L)
	u2=u.*(x.>=-L).*(x.<=L)+u₀(x.+L,2).*(x.<-L)+u₀(x.-L,2).*(x.>L)
	v2=u2-param.μ/3*Dx(h.^3 .*Dx(u2))./h

	η3=η.*(x.>=-L).*(x.<=L)+η₀(x.+L,3).*(x.<-L)+η₀(x.-L,3).*(x.>L)
	u3=u.*(x.>=-L).*(x.<=L)+u₀(x.+L,3).*(x.<-L)+u₀(x.-L,3).*(x.>L)
	v3=u3-param.μ/3*Dx(h.^3 .*Dx(u3))./h

	init1=Init(mesh,η1,v1)
	init2=Init(mesh,η2,v2)
	init3=Init(mesh,η3,v3)

	param2=merge(param,(L=-mesh.xmin,T = 8000, dt = 0.05, ns=1000))
	@time model=WhithamGreenNaghdi(param2;SGN=true,dealias=0,precond=false, gtol=1e-12)

	problem1 = Problem(model, init1, param2)
	problem2 = Problem(model, init2, param2)
	problem3 = Problem(model, init3, param2)
	solve!( problem1 )
	solve!( problem2 )
	solve!( problem3 )


	p = plot(layout=(2,1))
	plot_solution!(p,problem1,t=2000)
	plot_solution!(p,problem2,t=2000)
	plot_solution!(p,problem3,t=2000)
	savefig("ColdHot2000.pdf")
	xlims!(p[1,1],-L-200,-L+50)
	savefig("ColdHot2000Xzoom.pdf")
	ylims!(p[1,1],0.095,0.1)
	savefig("ColdHot2000Yzoom.pdf")

	p = plot(layout=(2,1))
	plot_solution!(p,problem1,t=4000)
	plot_solution!(p,problem2,t=4000)
	plot_solution!(p,problem3,t=4000)
	savefig("ColdHot4000.pdf")
	xlims!(p[1,1],-L-200,-L+50)
	savefig("ColdHot4000Xzoom.pdf")
	ylims!(p[1,1],0.095,0.1)
	savefig("ColdHot4000Yzoom.pdf")


	p = plot(layout=(2,1))
	plot_solution!(p,problem1,t=8000)
	plot_solution!(p,problem2,t=8000)
	plot_solution!(p,problem3,t=8000)
	savefig("ColdHot8000.pdf")
	xlims!(p[1,1],-L-200,-L+50)
	savefig("ColdHot8000Xzoom.pdf")
	ylims!(p[1,1],0.095,0.1)
	savefig("ColdHot8000Yzoom.pdf")


	display(p)

	save(problem1,"ColdHot1")
	save(problem2,"ColdHot2")
	save(problem3,"ColdHot3")
end

#---- Wavebreaking for the quasilinear Boussinesq-Whitham equation
function Wavebreaking1()
	param = ( μ  = 1/2,
			ϵ  = 1,
        	N  = 2^18,
            L  = 3,
            T  = .04,
            dt = 0.00001)
	η = x -> 2 .^(-abs.((x).^4))
	v = x -> -100*x.*η(x)
	init = Init(η,v)
	model = WhithamBoussinesq(merge(param,(α=1/2,)))
	problem = Problem(model, init, param)
	@time solve!( problem )
	p = plot(layout=(2,1))
	plot_solution!( p, problem, .03898)
	plot!(xlims=[(0.0018,0.0025) (-100000,100000 )])
	display(p)
end

#---- Wavebreaking for the quasilinear Boussinesq-Whitham equation ?
function Wavebreaking2()
	param = ( μ  = 1/2,
				ϵ  = 1,
	        	N  = 2^17, # you need to add more modes
	            L  = 3,
	            T  = 1,
	            dt = 0.0001) # you need to lower this
	η = x -> 2 .^(-abs.((x).^4))
	v = x -> -10*x.*η(x)
	init = Init(η,v)
	model = WhithamBoussinesq(merge(param,(α=1/2,)))
	problem = Problem(model, init, param)
	@time solve!( problem )
	p = plot(layout=(2,1))
	plot_solution!( p, problem,t= .92)
	plot!(xlims=[(0.275,0.285) (-100000,100000 )])
	display(p)
end
