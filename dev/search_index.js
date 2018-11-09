var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Documentation",
    "title": "Documentation",
    "category": "page",
    "text": ""
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
    "location": "basics/#Abstract-types-1",
    "page": "Code basics",
    "title": "Abstract types",
    "category": "section",
    "text": "TimeSolver (RK4, Euler, etc), \nAbstractModel (Cheng, Matsuno, etc), \nInitialData (Bump, SolitaryWave, etc)Instances are created from Parameters type.ParametersUne structure Problem  représente un problème donné que l\'on va résoudre. Les données seront stockées dans data, qui est vide initialement.Problem"
},

{
    "location": "basics/#Donnee-initiale-1",
    "page": "Code basics",
    "title": "Donnee initiale",
    "category": "section",
    "text": "mutable struct Bump <: InitialData  \n\n    # est-ce utile que ce soit mutable ?\n    # il faudrait forcer que toutes les InitialData (et pas seulement Bump) donnent h et u vecteurs réels\n\n    h   :: Array{Float64,1}\n    u   :: Array{Float64,1}\n\n    function Bump(p :: Parameters) # Une donnée initiale est construite à partir de son nom et de la donnée de Parameters\n    	mesh  = Mesh(-p.L, p.L, p.N)\n    	h = exp.(-(mesh.x).^2)\n    	u  = zeros(Complex{Float64}, mesh.N)\n    	new(h,u)\n    end\nend"
},

{
    "location": "basics/#Modele-1",
    "page": "Code basics",
    "title": "Modele",
    "category": "section",
    "text": "mutable struct Matsuno <: Model  \n\n    # Données qui seront utilisées dans les fonctions init, build, fwave! (et fig pour label)\n    label   :: String\n    h       :: Array{Complex{Float64},1}\n    u	      :: Array{Complex{Float64},1}\n    Gamma   :: Array{Float64,1}\n    Dx      :: Array{Complex{Float64},1}\n    H       :: Array{Complex{Float64},1}\n    Pi      :: BitArray{1}\n    Px      :: FFTW.cFFTWPlan{Complex{Float64},-1,false,1}\n    epsilon :: Float64\n    hnew    :: Array{Complex{Float64},1}\n    unew    :: Array{Complex{Float64},1}\n    Int1    :: Array{Complex{Float64},1}\n    Int2    :: Array{Complex{Float64},1}\n    Int3    :: Array{Complex{Float64},1}\n\n    function Matsuno(par::Parameters) # Un modèle est construit à partir de son nom et de la donnée de Parameters\n        label = \"Matsuno\"\n        mesh = Mesh(-par.L, par.L, par.N)\n        Gamma = abs.(mesh.k)\n        epsilon= par.epsilon\n        Dx    =  1im * mesh.k            # Differentiation\n        H     = -1im * sign.(mesh.k)     # Hilbert transform\n        Pi    = Gamma .< mesh.kmax * 2/3 # Dealiasing low-pass filter\n        h0 = zeros(Complex{Float64}, par.N)\n        Px  = plan_fft(h0; flags = FFTW.MEASURE)\n        h, u, hnew, unew ,Int1, Int2, Int3 = similar(h0), similar(h0), similar(h0), similar(h0), similar(h0), similar(h0), similar(h0)\n        new(label, h, u, Gamma, Dx, H, Pi, Px, epsilon, hnew, unew, Int1, Int2, Int3)\n    end\nend"
},

{
    "location": "basics/#Fonctions-1",
    "page": "Code basics",
    "title": "Fonctions",
    "category": "section",
    "text": "init et build pour le modèle Matsuno (c\'est le même pour Cheng)Pour définir complètement un nouveau modèle, il faut définir la structure au dessus, les fonctions init et build au dessous, ainsi que la fonction fwave (que je n\'ai pas recopiée)function init(m :: Matsuno, data::InitialData)\n\n         return (m.Pi.*fft(data.h),m.Pi.*fft(data.u))\nend\n\nfunction build(m :: Matsuno, h :: Array{Complex{Float64},1}, u :: Array{Complex{Float64},1})\n\n         return InitialData(real(ifft(h)),real(ifft(u)))\nend\n\n\nfunction solve!(p::Problem)\n\n    model = p.model( p.param )    # définit en particulier init et fwave! utilisés ci-dessous\n  	(h,u) = init(model,p.initial) # La fonction init, définie pour chaque modèle ransforme la donnée initiale\n  	solver = p.solver( p.param.N )\n\n  	times = Times(p.param.dt, p.param.T)\n  	prog = Progress(times.Nt,1) # progress bar\n\n    push!(p.data,(h,u))\n\n    for l in range(1,times.Nt-1)\n\n        step!( solver, model, fwave!, h, u, p.param.dt)   # la fonction step ne change pas par rapport aux versions précédentes\n        push!( p.data, (h,u))   \n\n        next!(prog)\n    end\n\nend"
},

{
    "location": "basics/#Document-principal-(ce-que-voit-l\'utilisateur)-1",
    "page": "Code basics",
    "title": "Document principal (ce que voit l\'utilisateur)",
    "category": "section",
    "text": "epsilon = 1/2\nN       = 2^12\nL       = 10\nT       = 5\ndt      = 0.001\n\nparam = Parameters(epsilon,N,L,T,dt)\n\nproblems = [Problem(Cheng,Bump,param,RK4),Problem(Matsuno,Bump,param,RK4)]Cheng est le modèle utilisé. Il prend en valeur Bump, param et définit 3 fonctioninit qui construit une donnée initiale à partir de Bump : Uinit=init(Bump,param)\nFwave utilisée pour résoudre dt U= Fwave(U) (avec donnée initiale Uinit)\nbuild qui reconstruit la donnée finale (c\'est l\'application inverse de init). Ufin=final(U,param)Bump est la donnée initiale, on pourrait stocker des exemples pertinents dans un dossier à part (comme pour les modèles et les solvers)param doit être défini par l\'utilisateur : c\'est grâce à lui que l\'on reconstruit mesh, time, etc. dans les différentes fonctions.Il pourrait changer entre deux simulation  :simuls = [(Cheng,Bump,param1,RK4),(Cheng,Bump,param2,RK4)]"
},

{
    "location": "basics/#Le-solver-RK4-1",
    "page": "Code basics",
    "title": "Le solver RK4",
    "category": "section",
    "text": "ne change pas par rapport aux versions précédentes. Il n\'a besoin que de param.N en argument (on pourrait par souci de cohérence lui donner carrément param en argument)for p in problems\n    solve!(p)\nendLe résultat de la simulation est stocké dans p.data pour chaque élément de problems il faut ensuite faire un fig pour plotter les solutions (c\'est là que la fonction build des modèles sera utile)"
},

]}
