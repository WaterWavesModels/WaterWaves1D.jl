using Test

#--- parameters
para = ( ϵ  = 0.1, μ = 0.1)  # physical parameters
paraX= ( N  = 6, L  = 3)   # mesh with 6 collocation points: x=[-3.0, -2.0, -1.0, 0.0, 1.0, 2.0]
paraT= ( T  = 2.5, dt = 1) # timegrid with three instants: t=[0.0, 1.0, 2.0]
param = merge(para,paraX)  # used to construct models
parap = merge(paraX,paraT) # used to construct problems

#--- initial data
init     = Init(x->exp.(-x.^2),x-> x )

#--- models
models=[]
push!(models,Boussinesq(param;
                    a=-1,b=3,
                    dealias=2,ktol=1e-12,
                    label="Boussinesq") )

push!(models,WhithamBoussinesq(param;
				α = 1.5, a = -1, b = 3,
				dealias = 2,
				ktol	= 1e-12,
				label 	= "Whitham-Boussinesq"
				) )

push!(models,WhithamBoussinesq(param;
				Boussinesq=true,
				α = 1.5, a = -1, b = 3,
				dealias = 2,
				ktol	= 1e-12,
				label 	= "(Whitham-)Boussinesq"
				) )

push!(models,DeepQuadratic_fast( param;
                    dealias=true, label="fast deep quadratic") )

push!(models,DeepQuadratic( param;
                    dealias=false, label="naive deep quadratic") )

push!(models,Matsuno_fast( param;
                    dealias=true, label="fast Matsuno") )

push!(models,Matsuno( param;
                    dealias=false, label="naive Matsuno") )


push!(models,IsobeKakinuma(param;
            dealias = 2,
            ktol	= 1e-10,
            iterate = true,
            gtol	= 1e-14,
            precond = false,
            restart	= 1,
            maxiter	= 2,
            label	= "Isobe-Kakinuma with GMRES"
			) )

push!(models, modifiedMatsuno(param;
			ν=2,ktol=01e-10,dealias=2,
			label="modified Matsuno") )

push!(models,IsobeKakinuma(param::NamedTuple;
			dealias = 0,
			ktol	= 0,
			iterate = false,
			gtol	= 0,
			precond = false,
			restart	= nothing,
			maxiter	= nothing,
			label	= "Isobe-Kakinuma with LU"
			) )

push!(models,NonHydrostatic(param;
            dealias = 2,
            ktol	= 1e-10,
            iterate = true,
            gtol	= 1e-14,
            precond = false,
            restart	= 1,
            maxiter	= 2,
            label	= "non-hydrostatic with GMRES"
			) )

push!(models,NonHydrostatic(param::NamedTuple;
			dealias = 0,
			ktol	= 0,
			iterate = false,
			gtol	= 0,
			precond = false,
			restart	= nothing,
			maxiter	= nothing,
			label	= "non-hydrostatic with LU"
			) )

push!(models,SquareRootDepth(param;
            dealias = 2,
            ktol	= 1e-10,
            iterate = true,
            gtol	= 1e-14,
            precond = false,
            restart	= 1,
            maxiter	= 2,
            label	= "square-root depth with GMRES"
			) )

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

push!(models,SerreGreenNaghdi(param;
            dealias = 2,
            ktol	= 1e-10,
            iterate = true,
            gtol	= 1e-14,
            precond = false,
            restart	= 1,
            maxiter	= 2,
            label	= "Green-Naghdi with GMRES"
			) )

push!(models,SerreGreenNaghdi(param::NamedTuple;
			dealias = 0,
			ktol	= 0,
			iterate = false,
			gtol	= 0,
			precond = false,
			restart	= nothing,
			maxiter	= nothing,
			label	= "Green-Naghdi with LU"
			) )

WhithamGreenNaghdi

push!(models,WhithamGreenNaghdi(param;
            dealias = 2,
            ktol	= 1e-10,
            iterate = true,
            gtol	= 1e-14,
            precond = false,
            restart	= 1,
            maxiter	= 2,
            label	= "Whitham-Green-Naghdi with GMRES"
			) )

push!(models,WhithamGreenNaghdi(param::NamedTuple;
			dealias = 0,
			ktol	= 0,
			iterate = false,
			gtol	= 0,
			precond = false,
			restart	= nothing,
			maxiter	= nothing,
			label	= "Whitham-Green-Naghdi with LU"
			) )

push!(models,WhithamGreenNaghdi(param;
			SGN=true,
            dealias = 2,
            ktol	= 1e-10,
            iterate = true,
            gtol	= 1e-14,
            precond = false,
            restart	= 1,
            maxiter	= 2,
            label	= "(Whitham-)Green-Naghdi with GMRES"
			) )

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

push!(models, WWn(param;
			ν		= 2,
			n		= 1,
			δ		= 0.1,
			m		= -12,
			ktol	= 1e-12,
			dealias	= 2,
			label	= "WW1"
			) )

push!(models, WWn(param;
			ν		= 1/2,
			n		= 7,
			δ		= 0.1,
			m		= -Inf,
			ktol	= 1e-12,
			dealias	= 2,
			label	= "WW2"
			) )

push!(models, WWn(param;
			ν		= 0,
			n		= 3,
			δ		= 0.1,
			m		= (-1,2),
			ktol	= 1e-12,
			dealias	= 2,
			label	= "WW3"
			) )

push!(models, WaterWaves(param;
					ν	    = 1,
					IL	    = true,
					method  = 1,
					tol	    = 1e-8,
					maxiter = 5,
					dealias	= 2,
					ktol	= 1e-10,
					label	= "water waves 1",
					verbose	= false) )


push!(models, WaterWaves(param;
					ν	    = 2,
					IL	    = false,
					method  = 2,
					tol	    = 1e-8,
					maxiter = 5,
					dealias	= 1,
					ktol	= 1e-10,
					label	= "water waves 2",
					verbose	= false) )


push!(models, WaterWaves(param;
					ν	    = 1/2,
					IL	    = false,
					method  = 3,
					tol	    = 1e-8,
					maxiter = 5,
					dealias	= 0,
					ktol	= 1e-10,
					label	= "water waves 3",
					verbose	= false) )


#--- tests
for model in models
    # build the initial-value problem
    problem = Problem( model, init, parap )
    # solve the initial-value problem
    solve!(problem ; verbose=false)
    # check everything went well
    @testset "Initial value problem with the model $(model.label)" begin
        @test length(problem.data.U)==length(Times(paraT).tc)
        @test size(problem.data.U[end])==(length(Mesh(paraX).x),2)
        @test !any(isnan,problem.data.U[end][1])
        @test !any(isnan,problem.data.U[end][2])
    end
end
