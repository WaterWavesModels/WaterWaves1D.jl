#md # # Animation
#md #
#md # deep water problem solved with Cheng model animation
#md #
#md # [`notebook`](@__NBVIEWER_ROOT_URL__notebooks/animation.ipynb)

#using DeepWaterModels
include("../src/includeall.jl")
using FFTW
using Plots
using ProgressMeter

#----

param = Parameters( ϵ  = 1/2,
                    N  = 2^12,
                    L  = 10,
                    T  = 5.0,
                    dt = 0.001)

initial = BellCurve(param,2.5)
solver  = RK4(param)
model   = CGBSW(param)
problem = Problem( model, initial, param, solver )

#----

function create_animation( p::Problem )

    h = construct(p.model,p.initial)[1]
    u = construct(p.model,p.initial)[2]

    prog = Progress(p.times.Nt,1)

    hr = real(similar(h))

    anim = @animate for l in range(1,p.times.Nt-1)

        dt = p.times.t[l+1]-p.times.t[l]

        step!(p.solver, p.model, h, u, dt)

        pl = plot(layout=(2,1))

        hr = real(ifft(h))

        plot!(pl[1,1], p.mesh.x, hr;
              ylims=(-0.6,1),
              title="physical space",
              label=p.model.label)

        plot!(pl[2,1], fftshift(p.mesh.k),
              log10.(1e-18.+abs.(fftshift(h)));
              title="frequency",
              label=p.model.label)

        next!(prog)

    end when mod(l, 200) == 0

    gif(anim, "anim.gif", fps=15); nothing

end

create_animation( problem )

#----
#md # ![](anim.gif)
