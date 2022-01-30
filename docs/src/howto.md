# How-to guides

## Add your model

## Add your initial data

## Add your solver






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

## Shallow water models

```@docs
Cheng
```

```@docs
Matsuno
```

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

Cheng est le modèle utilisé. Il prend en valeur Bump, param et définit 3 fonction:
  1. `init` qui construit une donnée initiale à partir de Bump : Uinit=init(Bump,param)
  2. `Fwave` utilisée pour résoudre dt U= Fwave(U) (avec donnée initiale Uinit)
  3. `build` qui reconstruit la donnée finale (c'est l'application inverse de init). Ufin=final(U,param)

```@docs
init(::Cheng)
init(::Matsuno)
```

```@docs
build(::Cheng)
build(::Matsuno)
```

```julia
simuls = [(Cheng,Bump,param1,RK4),(Cheng,Bump,param2,RK4)]
```

## RK4 solver

```@docs
RK4
```

```julia
for p in problems
    solve!(p)
end
```  
