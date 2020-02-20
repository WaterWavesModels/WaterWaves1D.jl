#md # # Animation
#md #
#md # shallow water problem solved with Boussinesq-Whitham model animation
#md #
#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/notebooks/animation.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/notebooks/animation.ipynb)

#using ShallowWaterModels
include("../src/dependencies.jl")

#----

param = ( μ  = 1/20,
          ϵ  = 1/2,
          N  = 2^11,
          L  = 10,
          T  = 8.0,
          dt = 0.001,
          α  = 1,
          θ  = 2)

initial = BellCurveExplicit(param)
solver  = RK4(param)
model   = fdBoussinesq(param)
problem = Problem( model, initial, param )

#----
print("\nNow solving the model ",problem.model.label,"\n")
@time solve!( problem )
#----

print("\nNow generating the animation\n")
@time create_animation( problem )

#----
#md # ![](anim.gif)
