var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Documentation",
    "title": "Documentation",
    "category": "page",
    "text": ""
},

{
    "location": "#DeepWaterModels.BellCurve",
    "page": "Documentation",
    "title": "DeepWaterModels.BellCurve",
    "category": "type",
    "text": "BellCurve(param,theta)\n\nh = 2^(-x^theta)\n\nu = 0\n\n\n\n\n\n"
},

{
    "location": "#DeepWaterModels.CGBSW",
    "page": "Documentation",
    "title": "DeepWaterModels.CGBSW",
    "category": "type",
    "text": "CGBSW( params )\n\n\n\n\n\n"
},

{
    "location": "#DeepWaterModels.CGBSW_naive",
    "page": "Documentation",
    "title": "DeepWaterModels.CGBSW_naive",
    "category": "type",
    "text": "CGBSW_naive( params )\n\n\n\n\n\n"
},

{
    "location": "#DeepWaterModels.HighFreq",
    "page": "Documentation",
    "title": "DeepWaterModels.HighFreq",
    "category": "type",
    "text": "HighFreq(param,s,k)\n\n\n\n\n\n"
},

{
    "location": "#DeepWaterModels.Matsuno_mod_naive",
    "page": "Documentation",
    "title": "DeepWaterModels.Matsuno_mod_naive",
    "category": "type",
    "text": "Modified Matsuno models with a naive step function\nMatsuno_mod_naive(params)\n\n\n\n\n\n"
},

{
    "location": "#DeepWaterModels.Matsuno_naive",
    "page": "Documentation",
    "title": "DeepWaterModels.Matsuno_naive",
    "category": "type",
    "text": "Matsuno(params)\n\n\n\n\n\n"
},

{
    "location": "#DeepWaterModels.Random",
    "page": "Documentation",
    "title": "DeepWaterModels.Random",
    "category": "type",
    "text": "Random(param,s,k)\n\n\n\n\n\n"
},

{
    "location": "#DeepWaterModels.mapfro-Tuple{CGBSW,Array{Complex{Float64},2}}",
    "page": "Documentation",
    "title": "DeepWaterModels.mapfro",
    "category": "method",
    "text": "mapfro(CGBSW, data)\n\n\n\n\n\n"
},

{
    "location": "#DeepWaterModels.mapfro-Tuple{CGBSW_naive,Array{Complex{Float64},2}}",
    "page": "Documentation",
    "title": "DeepWaterModels.mapfro",
    "category": "method",
    "text": "mapfro(CGBSW_naive, data)\n\n\n\n\n\n"
},

{
    "location": "#DeepWaterModels.mapfro-Tuple{Matsuno,Array{Complex{Float64},2}}",
    "page": "Documentation",
    "title": "DeepWaterModels.mapfro",
    "category": "method",
    "text": "mapfro(Matsuno, data)\n\n\n\n\n\n"
},

{
    "location": "#DeepWaterModels.mapfro-Tuple{Matsuno_mod_naive,Array{Complex{Float64},2}}",
    "page": "Documentation",
    "title": "DeepWaterModels.mapfro",
    "category": "method",
    "text": "mapfro(Matsuno, data)\n\n\n\n\n\n"
},

{
    "location": "#DeepWaterModels.mapfro-Tuple{Matsuno_naive,Array{Complex{Float64},2}}",
    "page": "Documentation",
    "title": "DeepWaterModels.mapfro",
    "category": "method",
    "text": "mapfro(Matsuno, data)\n\n\n\n\n\n"
},

{
    "location": "#DeepWaterModels.mapto-Tuple{CGBSW,InitialData}",
    "page": "Documentation",
    "title": "DeepWaterModels.mapto",
    "category": "method",
    "text": "mapto(CGBSW, data)\n\n\n\n\n\n"
},

{
    "location": "#DeepWaterModels.mapto-Tuple{CGBSW_naive,InitialData}",
    "page": "Documentation",
    "title": "DeepWaterModels.mapto",
    "category": "method",
    "text": "mapto(CGBSW_naive, data)\n\n\n\n\n\n"
},

{
    "location": "#DeepWaterModels.mapto-Tuple{Matsuno,InitialData}",
    "page": "Documentation",
    "title": "DeepWaterModels.mapto",
    "category": "method",
    "text": "mapto(Matsuno, data)\n\n\n\n\n\n"
},

{
    "location": "#DeepWaterModels.mapto-Tuple{Matsuno_mod_naive,InitialData}",
    "page": "Documentation",
    "title": "DeepWaterModels.mapto",
    "category": "method",
    "text": "mapto(Matsuno, data)\n\n\n\n\n\n"
},

{
    "location": "#DeepWaterModels.mapto-Tuple{Matsuno_naive,InitialData}",
    "page": "Documentation",
    "title": "DeepWaterModels.mapto",
    "category": "method",
    "text": "mapto(Matsuno, data)\n\n\n\n\n\n"
},

{
    "location": "#DeepWaterModels.jl-Documentation-1",
    "page": "Documentation",
    "title": "DeepWaterModels.jl Documentation",
    "category": "section",
    "text": "Modules = [DeepWaterModels]\nOrder   = [:type, :function]"
},

{
    "location": "basics/#",
    "page": "Code basics",
    "title": "Code basics",
    "category": "page",
    "text": ""
},

{
    "location": "basics/#Code-basics-1",
    "page": "Code basics",
    "title": "Code basics",
    "category": "section",
    "text": ""
},

{
    "location": "basics/#DeepWaterModels.Problem",
    "page": "Code basics",
    "title": "DeepWaterModels.Problem",
    "category": "type",
    "text": "Problem( model, initial, param, solver)\n\nmodel   : CGBSW or Matsuno\ninitial : BellCurve\nparam   : must contain N, L, T, dt for Mesh and Times, may contain additional data for Models (ϵ)\nsolver  : RK4 (optional)\n\n\n\n\n\n"
},

{
    "location": "basics/#Abstract-types-1",
    "page": "Code basics",
    "title": "Abstract types",
    "category": "section",
    "text": "TimeSolver (RK4, Euler, etc), \nAbstractModel (Cheng, Matsuno, etc), \nInitialData (Bump, SolitaryWave, etc)Instances are created from Parameters type.ParametersUne structure Problem  représente un problème donné que l\'on va résoudre. Les données seront stockées dans data, qui est vide initialement.Problem"
},

{
    "location": "basics/#Initial-data-1",
    "page": "Code basics",
    "title": "Initial data",
    "category": "section",
    "text": "Bump"
},

{
    "location": "basics/#DeepWaterModels.Matsuno",
    "page": "Code basics",
    "title": "DeepWaterModels.Matsuno",
    "category": "type",
    "text": "Matsuno(params)\n\n\n\n\n\n"
},

{
    "location": "basics/#Deep-water-models-1",
    "page": "Code basics",
    "title": "Deep water models",
    "category": "section",
    "text": "ChengMatsunosolve!(::Problem)"
},

{
    "location": "basics/#Main-program-1",
    "page": "Code basics",
    "title": "Main program",
    "category": "section",
    "text": "epsilon = 1/2\nN       = 2^12\nL       = 10\nT       = 5\ndt      = 0.001\n\nparam = Parameters(epsilon,N,L,T,dt)\n\nproblems = [Problem(Cheng,Bump,param,RK4),Problem(Matsuno,Bump,param,RK4)]Cheng est le modèle utilisé. Il prend en valeur Bump, param et définit 3 fonction:init qui construit une donnée initiale à partir de Bump : Uinit=init(Bump,param)\nFwave utilisée pour résoudre dt U= Fwave(U) (avec donnée initiale Uinit)\nbuild qui reconstruit la donnée finale (c\'est l\'application inverse de init). Ufin=final(U,param)init(::Cheng)\ninit(::Matsuno)build(::Cheng)\nbuild(::Matsuno)simuls = [(Cheng,Bump,param1,RK4),(Cheng,Bump,param2,RK4)]"
},

{
    "location": "basics/#DeepWaterModels.RK4",
    "page": "Code basics",
    "title": "DeepWaterModels.RK4",
    "category": "type",
    "text": "RK4(params)\n\nRunge-Kutta fourth order solver.\n\n\n\n\n\n"
},

{
    "location": "basics/#RK4-solver-1",
    "page": "Code basics",
    "title": "RK4 solver",
    "category": "section",
    "text": "RK4for p in problems\n    solve!(p)\nend"
},

{
    "location": "examples/animation/#",
    "page": "Animation",
    "title": "Animation",
    "category": "page",
    "text": "EditURL = \"https://github.com/WaterWavesModels/DeepWaterModels.jl/blob/master/examples/animation.jl\""
},

{
    "location": "examples/animation/#Animation-1",
    "page": "Animation",
    "title": "Animation",
    "category": "section",
    "text": "deep water problem solved with Cheng model animationnotebook#using DeepWaterModels\ninclude(\"../src/dependencies.jl\")param = ( ϵ  = 1/2,\n          N  = 2^12,\n          L  = 10,\n          T  = 5.0,\n          dt = 0.001)\n\ninitial = BellCurve(param,2.5)\nsolver  = RK4(param)\nmodel   = Matsuno(param)\nproblem = Problem( model, initial, param )print(\"\\nNow solving the model \",problem.model.label,\"\\n\")\n@time solve!( problem )print(\"\\nNow generating the animation\\n\")\n@time create_animation( problem )(Image: )This page was generated using Literate.jl."
},

{
    "location": "examples/two_problems/#",
    "page": "Example",
    "title": "Example",
    "category": "page",
    "text": "EditURL = \"https://github.com/WaterWavesModels/DeepWaterModels.jl/blob/master/examples/two_problems.jl\""
},

{
    "location": "examples/two_problems/#Two-deep-water-problems-1",
    "page": "Example",
    "title": "Two deep water problems",
    "category": "section",
    "text": "notebook#using DeepWaterModels\ninclude(\"../src/dependencies.jl\")param = ( ϵ  = 1/2,\n        	N  = 2^12,\n            L  = 10,\n            T  = 5,\n            dt = 0.001)\n\ninit     = BellCurve(param)\n\nmodel1    = CGBSW(param)\nsolver1   = RK4(param,model1)\nproblem1 = Problem(model1, init, param, solver1);\n\nmodel2  = Matsuno(param)\nsolver2   = RK4(param,model2)\nproblem2 = Problem(model2, init, param, solver2);p = plot(layout=(2,1))\n\nproblems = [ problem1, problem2 ]\n\nfor problem in problems\n	print(\"\\nNow solving the model \",problem.model.label,\"\\n\")\n   	@time solve!( problem )\n   	fig_problem!( p, problem )\n\nend\n\n#savefig(\"two_problems.png\"); nothing # hide\ndisplay(p)(Image: )This page was generated using Literate.jl."
},

{
    "location": "contents/#",
    "page": "Contents",
    "title": "Contents",
    "category": "page",
    "text": ""
},

{
    "location": "contents/#Contents-1",
    "page": "Contents",
    "title": "Contents",
    "category": "section",
    "text": ""
},

{
    "location": "contents/#Index-1",
    "page": "Contents",
    "title": "Index",
    "category": "section",
    "text": ""
},

]}
