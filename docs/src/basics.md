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

# Initial data

```@docs
Bump
```

## Deep water models

```@docs
Chen
```

```@docs
Matsuno
```

## Functions

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


```@docs
solve!(::Problem)
```

## Main program

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
