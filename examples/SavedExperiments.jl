# # Saved Experiments
#
#md # [`notebook`](@__NBVIEWER_ROOT_URL__notebooks/two_problems.ipynb)
#
using ShallowWaterModels,Plots,FFTW,Statistics;gr()
#include("../src/dependencies.jl")

#---- Cold-Hot states transition for Serre-Green-Naghdi, discussed in the paper of Sergey Gavrilyuk, Boniface Nkonga, Keh-Ming Shyue and Lev Truskinovsky
function ColdHot(scenario)
	if scenario !=1 && scenario !=2
		error("the scenario must be 1 or 2")
	end
	param = ( μ  = 1, ϵ  = 1, N  = 2^13,
			h₀=1.0962,h₁=1.1,h₂=1.2)
	@time (η,u,v,mesh,par)=CnoidalWaveWhithamGreenNaghdi(param;P=200,SGN=true)
	x=mesh.x
	h=1 .+ η
	m=-sqrt(param.h₀*param.h₁*param.h₂)
	if scenario == 1
		meanη=mean(h)-1
		meanv=m*mean(1 ./h)
	else
		meanη=1.09808-1
		meanv=m/1.09808
	end

	M=-m*mean(1 ./ h)/sqrt(mean(h)) #Froude number, close to 1
	L=80*par.λ

	v.-=mean(v)
	v.+=m*mean(1 ./h)

	interp1(x)=sech.(x/3).^2
	η1₀(x)=(meanη.+(0.2-meanη)*interp1(x))
	η1=η.*(x.>=-L).*(x.<=L)+η1₀.(x.+L).*(x.<-L)+η1₀.(x.-L).*(x.>L)
	v1₀(x)=(meanv .+(v[Int(end/2)]-meanv)*interp1(x))
	v1=v.*(x.>=-L).*(x.<=L)+v1₀(x.+L).*(x.<-L)+v1₀(x.-L).*(x.>L)
	#v1 .+=m*mean(1 ./h)


	interp2(x)=sech.(x/2).^2
	η2₀(x)=(meanη.+(0.2-meanη)*interp2(x))
	η2=η.*(x.>=-L).*(x.<=L)+η2₀.(x.+L).*(x.<-L)+η2₀.(x.-L).*(x.>L)
	v2₀(x)=(meanv .+(v[Int(end/2)]-meanv)*interp2(x))
	v2=v.*(x.>=-L).*(x.<=L)+v2₀(x.+L).*(x.<-L)+v2₀(x.-L).*(x.>L)
	#v2 .+=m*mean(1 ./h)


	interp3(x)=sech.(10*x).^2
	η3₀(x)=(meanη.+(0.2-meanη)*interp3(x))
	η3=η.*(x.>=-L).*(x.<=L)+η3₀.(x.+L).*(x.<-L)+η3₀.(x.-L).*(x.>L)
	v3₀(x)=(meanv .+(v[Int(end/2)]-meanv)*interp3(x))
	v3=v.*(x.>=-L).*(x.<=L)+v3₀(x.+L).*(x.<-L)+v3₀(x.-L).*(x.>L)
	#v3 .+=m*mean(1 ./h)

	init1=Init(mesh,η1,v1)
	init2=Init(mesh,η2,v2)
	init3=Init(mesh,η3,v3)
	param2=merge(param,(L=-mesh.xmin,T = 5000, dt = 0.05, ns=1000))
	@time model=WhithamGreenNaghdi(param2;SGN=true,dealias=1,precond=false, gtol=1e-12)
	problem1 = Problem(model, init1, param2)
	problem2 = Problem(model, init2, param2)
	problem3 = Problem(model, init3, param2)
	solve!( problem1 )
	solve!( problem2 )
	solve!( problem3 )

	p = plot(layout=(2,1))
	fig_problem!(p,problem1,1000)
	fig_problem!(p,problem2,1000)
	fig_problem!(p,problem3,1000)
	savefig("ColdHot1000.pdf")
	xlims!(p[1,1],-L-100,-L+50)
	savefig("ColdHot1000Xzoom.pdf")
	ylims!(p[1,1],0.0975,0.1)
	savefig("ColdHot1000Yzoom.pdf")


	p = plot(layout=(2,1))
	fig_problem!(p,problem1,3000)
	fig_problem!(p,problem2,3000)
	fig_problem!(p,problem3,3000)
	savefig("ColdHot3000.pdf")
	xlims!(p[1,1],-L-100,-L+50)
	savefig("ColdHot3000Xzoom.pdf")
	ylims!(p[1,1],0.0975,0.1)
	savefig("ColdHot3000Yzoom.pdf")

	p = plot(layout=(2,1))
	fig_problem!(p,problem1)
	fig_problem!(p,problem2)
	fig_problem!(p,problem3)
	savefig("ColdHot5000.pdf")
	xlims!(p[1,1],-L-100,-L+50)
	savefig("ColdHot5000Xzoom.pdf")
	ylims!(p[1,1],0.0975,0.1)
	savefig("ColdHot5000Yzoom.pdf")

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
	fig_problem!( p, problem, .03898)
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
	fig_problem!( p, problem, .92)
	plot!(xlims=[(0.275,0.285) (-100000,100000 )])
	display(p)
end
