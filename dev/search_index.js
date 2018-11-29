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
    "location": "#DeepWaterModels.construct-Tuple{CGBSW,InitialData}",
    "page": "Documentation",
    "title": "DeepWaterModels.construct",
    "category": "method",
    "text": "construct(CGBSW, data)\n\n\n\n\n\n"
},

{
    "location": "#DeepWaterModels.construct-Tuple{Matsuno,InitialData}",
    "page": "Documentation",
    "title": "DeepWaterModels.construct",
    "category": "method",
    "text": "construct( matsuno, data)\n\n\n\n\n\n"
},

{
    "location": "#DeepWaterModels.reconstruct-Tuple{CGBSW,Array{Complex{Float64},1},Array{Complex{Float64},1}}",
    "page": "Documentation",
    "title": "DeepWaterModels.reconstruct",
    "category": "method",
    "text": "reconstruct(CGBSW, h, u)\n\n\n\n\n\n"
},

{
    "location": "#DeepWaterModels.reconstruct-Tuple{Matsuno,Array{Complex{Float64},1},Array{Complex{Float64},1}}",
    "page": "Documentation",
    "title": "DeepWaterModels.reconstruct",
    "category": "method",
    "text": "reconstruct( matsuno, h, u)\n\n\n\n\n\n"
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
    "location": "basics/#DeepWaterModels.Parameters",
    "page": "Code basics",
    "title": "DeepWaterModels.Parameters",
    "category": "type",
    "text": "Parameters( ϵ  = 1/2,                 N  = 2^12,                L  = 10,                T  = 5,                dt = 0.001)\n\n\n\n\n\n"
},

{
    "location": "basics/#DeepWaterModels.Problem",
    "page": "Code basics",
    "title": "DeepWaterModels.Problem",
    "category": "type",
    "text": "Problem( model, initial, param, solver)\n\nmodel   : CGBSW or Matsuno\ninitial : BellCurve\nparam   : Mesh, Frequency, epsilon\nsolver  : RK4\n\n\n\n\n\n"
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
    "text": "deep water problem solved with Cheng model animationnotebookusing DeepWaterModels\nusing FFTW\nusing Plots\nusing ProgressMeterparam = Parameters( ϵ  = 1/2,\n                    N  = 2^12,\n                    L  = 10,\n                    T  = 5.0,\n                    dt = 0.001)\n\ninitial = BellCurve(param,2.5)\nsolver  = RK4(param)\nmodel   = CGBSW(param)\nproblem = Problem( model, initial, param, solver )function create_animation( p::Problem )\n\n    h = construct(p.model,p.initial)[1]\n    u = construct(p.model,p.initial)[2]\n\n    prog = Progress(p.times.Nt,1)\n\n    hr = real(similar(h))\n\n    anim = @animate for l in range(1,p.times.Nt-1)\n\n        dt = p.times.t[l+1]-p.times.t[l]\n\n        step!(p.solver, p.model, h, u, dt)\n\n        pl = plot(layout=(2,1))\n\n        hr = real(ifft(h))\n\n        plot!(pl[1,1], p.mesh.x, hr;\n              ylims=(-0.6,1),\n              title=\"physical space\",\n              label=p.model.label)\n\n        plot!(pl[2,1], fftshift(p.mesh.k),\n              log10.(1e-18.+abs.(fftshift(h)));\n              title=\"frequency\",\n              label=p.model.label)\n\n        next!(prog)\n\n    end when mod(l, 200) == 0\n\n    gif(anim, \"anim.gif\", fps=15); nothing\n\nend\n\ncreate_animation( problem )(Image: )This page was generated using Literate.jl."
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
    "text": "notebookusing DeepWaterModels\nusing FFTW\nusing PlotsPlot results functionfunction fig_problem!( p, problem::Problem )\n\n    s = 0\n    (hhat,uhat) = problem.data[end]\n    (hr,ur)     = (real(ifft((problem.model.Gamma.^s).*hhat)),\n		   real(ifft(uhat)))\n\n    plot!(p[1,1], problem.model.mesh.x, hr;\n		  title=\"physical space\",\n	          label=problem.model.label)\n\n    plot!(p[2,1], fftshift(problem.model.mesh.k),\n                  log10.(1e-18.+abs.(fftshift(hhat)));\n		  title=\"frequency\",\n    	          label=problem.model.label)\n\nendparam = Parameters( ϵ  = 1/2,\n                    N  = 2^12,\n                    L  = 10,\n                    T  = 5,\n                    dt = 0.001)\n\ninit     = BellCurve(param)\nsolver   = RK4(param)\n\nmodel    = CGBSW(param)\nproblem1 = Problem(model, init, param, solver)\n\nmatsuno  = Matsuno(param)\nproblem2 = Problem(matsuno, init, param, solver);p = plot(layout=(2,1))\n\nproblems = [ problem1, problem2 ]\n\nfor problem in problems\n\n   solve!( problem )\n   fig_problem!( p, problem )\n\nend\n\nsavefig(\"two_problems.png\"); nothing # hide(Image: )This page was generated using Literate.jl."
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
