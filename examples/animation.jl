#md # # Animation
#md #
#md # deep water problem solved with Cheng model animation
#md #
#md # [`notebook`](@__NBVIEWER_ROOT_URL__notebooks/animation.ipynb)

#using DeepWaterModels
include("../src/dependencies.jl")

#----

param = ( ϵ  = 1/2,
          N  = 2^12,
          L  = 10,
          T  = 5.0,
          dt = 0.001)

initial = BellCurve(param,2.5)
#solver  = RK4(param)
model   = CGBSW(param)
problem = Problem( model, initial, param )

#----
print("\nNow solving the model ",problem.model.label,"\n")
@time solve!( problem )
#----

print("\nNow generating the animation\n")
@time create_animation( problem )

#----
#md # ![](anim.gif)
