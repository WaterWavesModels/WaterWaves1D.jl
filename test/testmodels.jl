using Test


#--- initial data (smooth, support compact). Zero velocity for compatibility between models.
init     = Init(x->exp.(-1/8*x.^4),x-> 0. * x )

#--- Tests on shallow water models

para = ( ϵ  = 0.01, μ = 0.01)  # physical parameters
paraX= ( N  = 2^6, L  = 4)   # mesh with 64 collocation points on [-4,4]
paraT= ( T  = 1e-1, dt = 1e-2) # timegrid with 10 instants: t=[0.0:1.0:10.0]/100
param = merge(para,paraX)  # used to construct models

models=AbstractModel[];precisions=Real[]

# Build shallow water models

push!(models,SaintVenant(param;
                    dealias=0,ktol=1e-12,
                    label="Saint-Venant"
					) )
push!(precisions,para.μ)

push!(models,SaintVenant_fast(param;
                    dealias=0,ktol=0,
                    label="Saint-Venant"
					) )
push!(precisions,para.μ)

push!(models,SerreGreenNaghdi(param;
            dealias = 0,
            ktol	= 1e-10,
            iterate = true,
            gtol	= 1e-14,
            precond = true,
            restart	= 15,
            maxiter	= 25,
            label	= "Green-Naghdi with GMRES"
			) )
push!(precisions,para.μ^2)

push!(models,SerreGreenNaghdi(param::NamedTuple;
			dealias = 1,
			ktol	= 0,
			iterate = false,
			gtol	= 0,
			precond = false,
			restart	= nothing,
			maxiter	= nothing,
			label	= "Green-Naghdi with LU"
			) )
push!(precisions,para.μ^2)

push!(models,WhithamGreenNaghdi(param::NamedTuple;
			SGN=true,
			dealias = 0,
			ktol	= 0,
			iterate = false,
			gtol	= 0,
			precond = false,
			restart	= nothing,
			maxiter	= nothing,
			label	= "(Whitham-)Green-Naghdi with LU"
			) )
push!(precisions,para.μ^2)


push!(models,WhithamGreenNaghdi(param;
			SGN=true,
            dealias = false,
            ktol	= 1e-10,
            iterate = true,
            gtol	= 1e-14,
            precond = false,
            restart	= 15,
            maxiter	= 25,
            label	= "(Whitham-)Green-Naghdi with GMRES"
			) )
push!(precisions,para.μ^2)


push!(models,WhithamGreenNaghdi(param;
            dealias = 0,
            ktol	= 1e-10,
            iterate = true,
            gtol	= 1e-14,
            precond = false,
            restart	= 5,
            maxiter	= 15,
            label	= "Whitham-Green-Naghdi with GMRES") )
push!(precisions,para.μ^2*param.ϵ)


push!(models,WhithamGreenNaghdi(param::NamedTuple;
			dealias = 0,
			ktol	= 0,
			iterate = false,
			gtol	= 0,
			precond = false,
			restart	= nothing,
			maxiter	= nothing,
			label	= "Whitham-Green-Naghdi with LU") )
push!(precisions,para.μ^2*param.ϵ)


push!(models,IsobeKakinuma(param::NamedTuple;
			dealias = 0,
			ktol	= 0,
			iterate = false,
			gtol	= 0,
			precond = false,
			restart	= nothing,
			maxiter	= nothing,
			label	= "Isobe-Kakinuma with LU") )
push!(precisions,para.μ^3)

push!(models,IsobeKakinuma(param;
            dealias = 0,
            ktol	= 1e-14,
            iterate = true,
            gtol	= 1e-14,
            precond = true,
            restart	= 10,
            maxiter	= 20,
            label	= "Isobe-Kakinuma with GMRES") )
push!(precisions,para.μ^3)


push!(models,NonHydrostatic(param;
            dealias = 2,
            ktol	= 1e-14,
            iterate = true,
            gtol	= 1e-14,
            precond = false,
            restart	= 1,
            maxiter	= 2,
            label	= "non-hydrostatic with GMRES") )
push!(precisions,para.μ/5)

push!(models,NonHydrostatic(param::NamedTuple;
			dealias = 0,
			ktol	= 0,
			iterate = false,
			gtol	= 0,
			precond = false,
			restart	= nothing,
			maxiter	= nothing,
			label	= "non-hydrostatic with LU") )
push!(precisions,para.μ/5)

push!(models,SquareRootDepth(param;
            dealias = 2,
            ktol	= 1e-10,
            iterate = true,
            gtol	= 1e-14,
            precond = false,
            restart	= 15,
            maxiter	= 20,
            label	= "square-root depth with GMRES") )
push!(precisions,max(para.μ*para.ϵ,para.μ^2))


push!(models,SquareRootDepth(param::NamedTuple;
			dealias = 0,
			ktol	= 0,
			iterate = false,
			gtol	= 0,
			precond = false,
			restart	= nothing,
			maxiter	= nothing,
			label	= "square-root depth with LU"
			) )
push!(precisions,max(para.μ*para.ϵ,para.μ^2))


push!(models,Boussinesq(param;
                    a=-1//2,b=5//12,
                    dealias=0,ktol=1e-12,
                    label="Boussinesq"
					) )
push!(precisions,max(para.μ*para.ϵ/4,para.μ^2))


push!(models,WhithamBoussinesq(param;
				Boussinesq=true,
				α = 1, a = -1/2, b = 5/12,
				dealias = 0,
				ktol	= 1e-12,
				label 	= "(Whitham-)Boussinesq"
				) )
push!(precisions,max(para.μ*para.ϵ/4,para.μ^2))

push!(models,WhithamBoussinesq(param;
				α = 1.5, a = -1, b = 3,
				dealias = 0,
				ktol	= 1e-12,
				label 	= "Whitham-Boussinesq"
				) )
push!(precisions,para.μ*para.ϵ/4)

push!(models,Choi(param::NamedTuple;
			M=0,reg=true,
			dealias = 0,
			ktol	= 0,
			iterate = false,
			gtol	= 0,
			precond = false,
			restart	= nothing,
			maxiter	= nothing,
			label	= "Choi with LU and M=0") )
push!(precisions,para.μ)

push!(models,Choi(param;
			M=1,
            dealias = 1,
            ktol	= 1e-14,
            iterate = true,
            gtol	= 1e-14,
            precond = true,
            restart	= 10,
            maxiter	= 20,
            label	= "Choi with GMRES and M=1") )
push!(precisions,para.μ^2)

push!(models,Choi(param::NamedTuple;
			M=2,
			reg=true,
			dealias = 1,
			ktol	= 0,
			iterate = false,
			gtol	= 0,
			precond = false,
			restart	= nothing,
			maxiter	= nothing,
			label	= "Choi with LU and M=2") )
push!(precisions,para.μ^3)

paramrelax=merge(param,(a=100,))
push!(models,relaxedGreenNaghdi(paramrelax::NamedTuple;
			id=2,
			FG=true,
			dealias = 1,
			ktol	= 0,
			iterate = false,
			gtol	= 0,
			precond = false,
			restart	= nothing,
			maxiter	= nothing,
			label	= "relaxed Green-Naghdi (Favrie-Gavrilyuk)") )
push!(precisions,para.μ^2)

push!(models,relaxedGreenNaghdi(paramrelax::NamedTuple;
			id=1,
			FG=false,
			dealias = 1,
			ktol	= 1e-12,
            iterate = true,
            gtol	= 1e-14,
            precond = true,
            restart	= 10,
            maxiter	= 20,
			label	= "relaxed Green-Naghdi (Escalante-Dumbser-Castro)") )
push!(precisions,para.μ^2)

push!(models,Whitham(param::NamedTuple;
			KdV=false,
			BBM=false,
			dealias = 1,
			ktol	= 1e-12,
            label	= "Whitham") )
push!(precisions,para.μ*para.ϵ+para.ϵ^2)

push!(models,Whitham(param::NamedTuple;
			KdV=true,
			BBM=false,
			dealias = 1,
			ktol	= 1e-12,
            label	= "KdV") )
push!(precisions,para.μ^2+para.ϵ^2)

push!(models,Whitham(param::NamedTuple;
			KdV=false,
			BBM=true,
			dealias = 1,
			ktol	= 1e-12,
            label	= "BBM") )
push!(precisions,para.μ^2+para.ϵ^2)

push!(models,KdV(param::NamedTuple;
			dealias = 0,
			ktol	= 1e-12,
            label	= "KdV") )
push!(precisions,para.μ^2+para.ϵ^2)

push!(models,BBM(param::NamedTuple;
			dealias = 0,
			ktol	= 1e-12,
            label	= "BBM") )
push!(precisions,para.μ^2+para.ϵ^2)

# Build reference problem (water waves)
modelWW =  WaterWaves(merge(param,(ν=1,));
					IL	    = false,
					method  = 3,
					tol	    = 1e-8,
					maxiter = 5,
					dealias	= 0,
					ktol	= 1e-10,
					label	= "water waves (shallow water)",
					verbose	= false)
pbWW = Problem( modelWW, init, paraT )
solve!(pbWW;verbose=false)
ηWW,vWW,xWW=solution(pbWW)

for i in eachindex(models)
    # build the initial-value problem
    problem = Problem( models[i], init, paraT )
    # solve the initial-value problem
    solve!(problem ; verbose=false)
    # check the order of magnitude of the precision is correct (from above and from below)
    @testset "Initial value problem with the model $(models[i].label)" begin
		@test isapprox(solution(problem,x=xWW)[1] , ηWW, rtol = precisions[i]*paraT.T )
		@test !isapprox(solution(problem,x=xWW)[1] , ηWW, rtol = precisions[i]*paraT.T/10 )
    end
end

# Check SaintVenant is identical to SaintVenant_fast
problem1 = Problem( models[1], init, paraT )
problem2 = Problem( models[2], init, paraT )
# solve the initial-value problem
solve!(problem1 ; verbose=false)
solve!(problem2 ; verbose=false)
@test all(isapprox.(solution(problem1) , solution(problem2) ))

# Check two Boussinesq models are the same
model1=Boussinesq(param;
                    a=-1//2,b=5//12,
                    dealias=1,ktol=1e-12,
                    label="Boussinesq"
					) 
model2=WhithamBoussinesq(param;
				Boussinesq=true,
				a = -1/2, b = 5/12,
				dealias = 1,
				ktol	= 1e-12,
				label 	= "(Whitham-)Boussinesq"
				) 
problem1 = Problem( model1, init, paraT )
problem2 = Problem( model2, init, paraT )
# solve the initial-value problem
solve!(problem1 ; verbose=false)
solve!(problem2 ; verbose=false)
@test all(isapprox.(solution(problem1) , solution(problem2) ))


#--- Tests on deep water models
models=AbstractModel[];precisions=Real[]

para = ( ϵ  = 0.1, μ = 100)  # physical parameters
paraX= ( N  = 2^6, L  = 4)   # mesh with 64 collocation points on [-4,4]
paraT= ( T  = 1e-1, dt = 1e-2) # timegrid with 10 instants: t=[0.0:1.0:10.0]/100
param = merge(para,paraX)  # used to construct models

# Build deep layer models
push!(models, modifiedMatsuno(merge(param,(ν=1/√para.μ,));
			ktol=01e-14,dealias=0,
			label="modified Matsuno"
			) )
push!(precisions, para.ϵ.^2*para.μ)

push!(models,Airy(param;
                    label="Airy"
					) )
push!(precisions,para.ϵ*√para.μ)


push!(models, WWn(param;
			n		= 1,
			δ		= 0.1,
			m		= (-1,2),
			ktol	= 1e-12,
			dealias	= 1,
			label	= "WW1"
			) )
push!(precisions, para.ϵ*√para.μ)

push!(models, WWn(merge(param,(ν=1/√para.μ,));
			n		= 7,
			δ		= 0.1,
			m		= -Inf,
			ktol	= 1e-12,
			dealias	= 0,
			label	= "WW2"
			) )
push!(precisions, para.ϵ.^2*para.μ)

push!(models, WWn(param;
			n		= 3,
			δ		= 0.01,
			ktol	= 1e-14,
			dealias	= 0,
			label	= "WW3"
			) )
push!(precisions, para.ϵ.^3*para.μ^(3/2))

# Build reference problem (water waves)
modelWW = WaterWaves(param;
					IL	    = false,
					method  = 2,
					tol	    = 1e-12,
					maxiter = 25,
					dealias	= 0,
					ktol	= 1e-10,
					label	= "water waves (deep water)",
					verbose	= false)
pbWW = Problem( modelWW, init, paraT )
solve!(pbWW;verbose=false)
ηWW,vWW,xWW=solution(pbWW)


for i in eachindex(models)
    # build the initial-value problem
    problem = Problem( models[i], init, paraT )
    # solve the initial-value problem
    solve!(problem ; verbose=false)
    # check the order of magnitude of the precision is correct (from above and from below)
    @testset "Initial value problem with the model $(models[i].label)" begin
		@test isapprox(solution(problem,x=xWW)[1] , ηWW, rtol = precisions[i]*paraT.T/50 )
		@test !isapprox(solution(problem,x=xWW)[1] , ηWW, rtol = precisions[i]*paraT.T/500 )

    end
end

#--- Tests on infinite layer models
models=AbstractModel[];precisions=Real[]

para = ( ϵ  = 0.1, μ = Inf)  # physical parameters
paraX= ( N  = 2^6, L  = 4)   # mesh with 64 collocation points on [-4,4]
paraT= ( T  = 1e-1, dt = 1e-2) # timegrid with 10 instants: t=[0.0:1.0:10.0]/100
param = merge(para,paraX)  # used to construct models

# Build infinite layer models
push!(models, modifiedMatsuno(param;
			ktol=01e-14,dealias=0,
			label="modified Matsuno"
			) )
push!(precisions, para.ϵ.^2)

push!(models, WWn(param;
			n		= 1,
			δ		= 0.1,
			m		= (-1,2),
			ktol	= 1e-12,
			dealias	= 1,
			label	= "WW1") )
push!(precisions, para.ϵ)

push!(models, WWn(merge(param,(ν=0,));
			n		= 7,
			δ		= 0.1,
			m		= -Inf,
			ktol	= 1e-12,
			dealias	= 0,
			label	= "WW2") )
push!(precisions, para.ϵ.^2)


push!(models, WWn(param;
			IL		= true,
			n		= 3,
			δ		= 0.01,
			ktol	= 0,
			dealias	= 0,
			label	= "WW3") )
push!(precisions, para.ϵ.^3)



push!(models,AkersNicholls_fast( param;
                    dealias=true, label="fast deep quadratic") )
push!(precisions, para.ϵ.^2)

push!(models,AkersNicholls( param;
                    dealias=false, label="naive deep quadratic") )
push!(precisions, para.ϵ.^2)

push!(models,Matsuno_fast( param;
                    dealias=true, label="fast Matsuno") )
push!(precisions, para.ϵ.^2)

push!(models,Matsuno( param;
                    dealias=false, label="naive Matsuno") )
push!(precisions, para.ϵ.^2)


# Build reference problem (water waves)
modelWW =  WaterWaves(param;
					IL	    = true,
					method  = 1,
					tol	    = 1e-12,
					maxiter = 5,
					dealias	= 0,
					ktol	= 1e-14,
					label	= "water waves (infinite layer)",
					verbose	= false)
pbWW = Problem( modelWW, init, paraT )
solve!(pbWW;verbose=false)
ηWW,vWW,xWW=solution(pbWW)

for i in eachindex(models)
    # build the initial-value problem
    problem = Problem( models[i], init, paraT )
    # solve the initial-value problem
    solve!(problem ; verbose=false)
    # check the order of magnitude of the precision is correct (from above and from below)
    @testset "Initial value problem with the model $(models[i].label)" begin
		@test isapprox(solution(problem,x=xWW)[1] , ηWW, rtol = precisions[i]*paraT.T/50 )
		@test !isapprox(solution(problem,x=xWW)[1] , ηWW, rtol = precisions[i]*paraT.T/500 )

    end
end



#--- Tests on 2D models
@testset "2D models reduce to 1D models" begin

	# Initial data
	ζ(x,y) = 0.5 * cos.(π*x).*ones(length(y))';
	ux(x,y) = 0.5 * sin.(π*x).*ones(length(y))';
	uy(x,y) = zero(x).*ones(length(y))';

	init1D = Init(x->ζ(x,0),x-> ux(x,0));
	init2D = Init2D(ζ, ux, uy);

	# Build 1D problem
	model1D = SaintVenant(param; dealias = 0, smooth = 0)
	solver1D=RK4(model1D.mapto(init1D))
	problem1D = Problem(  model1D, init1D, paraT ; solver=solver1D )

	# Build 2D problem
	model2D=SaintVenant2D(param; dealias = 0 , hamiltonian = false, smooth=0)
	solver2D=RK4(model2D.mapto(init2D))
	problem2D = Problem(  model2D, init2D, paraT ; solver=solver2D )

	#---- Solve problems
	solve!(problem1D;verbose=false)
	solve!(problem2D;verbose=false)

	η,v= solution(problem1D);
	η2D,vx,vy,x,y = solution(problem2D);

	#---- Perform checks
	@test η2D[:,1] ≈ η
	@test sum(η2D,dims=2) ≈ param.N* η
	@test vx[:,1] ≈ v
	@test sum(vx,dims=2) ≈ param.N* v
	@test vy==zero(vy)
end

@testset "2D models are identical for irrotational data" begin

	# Initial data
	ζ(x,y) = 0.5 .*cos.(π*x).*cos.(π*y)';
	ux(x,y) = 0.5 .* cos.(π*y').*sin.(π*x);
	uy(x,y) = 0.5 .* sin.(π*y').*cos.(π*x);

	init = Init2D(ζ, ux, uy);

	# Build problems
	model=SaintVenant2D_fast(param; dealias = 1/4 , hamiltonian = false, smooth=1)
	model_hamiltonian=SaintVenant2D_fast(param; dealias = 1/4 , hamiltonian = true, smooth=1)

	solver=RK4(model.mapto(init))

	problem = Problem(  model, init, paraT ; solver=solver )
	problem_hamiltonian = Problem(  model_hamiltonian, init, paraT ; solver=solver )

	# Solve problems
	solve!(problem;verbose=false)
	solve!(problem_hamiltonian;verbose=false)

	# Perform checks
	η,vx,vy,x,y = solution(problem;T=1)
	hη,hvx,hvy,hx,hy = solution(problem_hamiltonian;T=1)

	@test hη ≈ η
	@test hvx ≈ vx
	@test hvy ≈ vy
end

@testset "2D models are not identical for non-irrotational data" begin

	# Initial data
	ζ(x,y) = 0.5 .*cos.(π*x).*cos.(π*y)';
	ux(x,y) = 0.5 .* cos.(π*y').*sin.(π*x);
	uy(x,y) = -0.5 .* sin.(π*y').*cos.(π*x);

	init = Init2D(ζ, ux, uy);

	# Build problems
	model=SaintVenant2D(param; dealias = 1/4 , hamiltonian = false, smooth=1)
	model_hamiltonian=SaintVenant2D(param; dealias = 1/4 , hamiltonian = true, smooth=1)

	solver=RK4(model.mapto(init))

	problem = Problem(  model, init, paraT ; solver=solver )
	problem_hamiltonian = Problem(  model_hamiltonian, init, paraT ; solver=solver )

	# Solve problems
	solve!(problem;verbose=false)
	solve!(problem_hamiltonian;verbose=false)

	# Perform checks
	η,vx,vy,x,y = solution(problem;T=1)
	hη,hvx,hvy,hx,hy = solution(problem_hamiltonian;T=1)

	@test !(hη ≈ η)
	@test !(hvx ≈ vx)
	@test !(hvy ≈ vy)
end