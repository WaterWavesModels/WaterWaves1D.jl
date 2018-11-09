# Code basics

## Abstract types

  - `TimeSolver` (`RK4`, `Euler`, etc), 
  - `AbstractModel` (`Cheng`, `Matsuno`, etc), 
  - `InitialData` (`Bump`, `SolitaryWave`, etc)

Instances are created from `Parameters` type.

```@docs
Parameters
```

Une structure `Problem`  représente un problème donné que l'on va résoudre.
Les données seront stockées dans data, qui est vide initialement.

```@docs
Problem
```

# Donnee initiale

```julia
mutable struct Bump <: InitialData  

    # est-ce utile que ce soit mutable ?
    # il faudrait forcer que toutes les InitialData (et pas seulement Bump) donnent h et u vecteurs réels

    h   :: Array{Float64,1}
    u   :: Array{Float64,1}

    function Bump(p :: Parameters) # Une donnée initiale est construite à partir de son nom et de la donnée de Parameters
    	mesh  = Mesh(-p.L, p.L, p.N)
    	h = exp.(-(mesh.x).^2)
    	u  = zeros(Complex{Float64}, mesh.N)
    	new(h,u)
    end
end
```

## Modele

```julia
mutable struct Matsuno <: Model  

    # Données qui seront utilisées dans les fonctions init, build, fwave! (et fig pour label)
    label   :: String
    h       :: Array{Complex{Float64},1}
    u	      :: Array{Complex{Float64},1}
    Gamma   :: Array{Float64,1}
    Dx      :: Array{Complex{Float64},1}
    H       :: Array{Complex{Float64},1}
    Pi      :: BitArray{1}
    Px      :: FFTW.cFFTWPlan{Complex{Float64},-1,false,1}
    epsilon :: Float64
    hnew    :: Array{Complex{Float64},1}
    unew    :: Array{Complex{Float64},1}
    Int1    :: Array{Complex{Float64},1}
    Int2    :: Array{Complex{Float64},1}
    Int3    :: Array{Complex{Float64},1}

    function Matsuno(par::Parameters) # Un modèle est construit à partir de son nom et de la donnée de Parameters
        label = "Matsuno"
        mesh = Mesh(-par.L, par.L, par.N)
        Gamma = abs.(mesh.k)
        epsilon= par.epsilon
        Dx    =  1im * mesh.k            # Differentiation
        H     = -1im * sign.(mesh.k)     # Hilbert transform
        Pi    = Gamma .< mesh.kmax * 2/3 # Dealiasing low-pass filter
        h0 = zeros(Complex{Float64}, par.N)
        Px  = plan_fft(h0; flags = FFTW.MEASURE)
        h, u, hnew, unew ,Int1, Int2, Int3 = similar(h0), similar(h0), similar(h0), similar(h0), similar(h0), similar(h0), similar(h0)
        new(label, h, u, Gamma, Dx, H, Pi, Px, epsilon, hnew, unew, Int1, Int2, Int3)
    end
end
```

## Fonctions

`init` et `build` pour le modèle Matsuno (c'est le même pour Cheng)

Pour définir complètement un nouveau modèle, il faut définir la structure au
dessus, les fonctions init et build au dessous, ainsi que la fonction `fwave`
(que je n'ai pas recopiée)

```julia
function init(m :: Matsuno, data::InitialData)

         return (m.Pi.*fft(data.h),m.Pi.*fft(data.u))
end

function build(m :: Matsuno, h :: Array{Complex{Float64},1}, u :: Array{Complex{Float64},1})

         return InitialData(real(ifft(h)),real(ifft(u)))
end


function solve!(p::Problem)

    model = p.model( p.param )    # définit en particulier init et fwave! utilisés ci-dessous
  	(h,u) = init(model,p.initial) # La fonction init, définie pour chaque modèle ransforme la donnée initiale
  	solver = p.solver( p.param.N )

  	times = Times(p.param.dt, p.param.T)
  	prog = Progress(times.Nt,1) # progress bar

    push!(p.data,(h,u))

    for l in range(1,times.Nt-1)

        step!( solver, model, fwave!, h, u, p.param.dt)   # la fonction step ne change pas par rapport aux versions précédentes
        push!( p.data, (h,u))   

        next!(prog)
    end

end
```

## Document principal (ce que voit l'utilisateur)

```julia
epsilon = 1/2
N       = 2^12
L       = 10
T       = 5
dt      = 0.001

param = Parameters(epsilon,N,L,T,dt)

problems = [Problem(Cheng,Bump,param,RK4),Problem(Matsuno,Bump,param,RK4)]
```

Cheng est le modèle utilisé. Il prend en valeur Bump, param et définit 3 fonction

  1. `init` qui construit une donnée initiale à partir de Bump : Uinit=init(Bump,param)
  2. `Fwave` utilisée pour résoudre dt U= Fwave(U) (avec donnée initiale Uinit)
  3. `build` qui reconstruit la donnée finale (c'est l'application inverse de init). Ufin=final(U,param)

Bump est la donnée initiale, on pourrait stocker des exemples pertinents dans un dossier à part (comme pour les modèles et les solvers)

param doit être défini par l'utilisateur : c'est grâce à lui que l'on reconstruit mesh, time, etc. dans les différentes fonctions.

Il pourrait changer entre deux simulation  :

```julia
simuls = [(Cheng,Bump,param1,RK4),(Cheng,Bump,param2,RK4)]
```

## Le solver RK4
ne change pas par rapport aux versions précédentes. Il n'a besoin que de param.N en argument (on pourrait par souci de cohérence lui donner carrément param en argument)

```julia
for p in problems
    solve!(p)
end
```  

Le résultat de la simulation est stocké dans p.data pour chaque élément de problems
il faut ensuite faire un fig pour plotter les solutions (c'est là que la fonction build des modèles sera utile)
