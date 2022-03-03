using Test


#--- initial data (smooth, support compact). Zero velocity for compatibility between models.
init     = Init(x->exp.(-1/8*x.^4),x-> 0. * x )

#--- Tests on shallow water models

para = ( ϵ  = 0.01, μ = 0.01)  # physical parameters
paraX= ( N  = 2^6, L  = 4)   # mesh with 64 collocation points on [-4,4]
paraT= ( T  = 1e-1, dt = 1e-2) # timegrid with 10 instants: t=[0.0:1.0:10.0]/100
param = merge(para,paraX)  # used to construct models

models=[];precisions=[]

# Build shallow water models
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

push!(models,SaintVenant(param;
                    dealias=0,ktol=1e-12,
                    label="Saint-Venant"
					) )
push!(precisions,para.μ)


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
				ktol	= 1e-14,
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

realp=[]
for i in 1:length(models)
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


#--- Tests on deep water models
models=[];precisions=[]

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


for i in 1:length(models)
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
models=[];precisions=[]

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

for i in 1:length(models)
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
