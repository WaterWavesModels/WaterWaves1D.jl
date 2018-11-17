#md # # Animation
#md
#md # deep water problem solved with Cheng model animation
#
#md # [`notebook`](@__NBVIEWER_ROOT_URL__notebooks/animation.ipynb)
#
using DeepWaterModels
using FFTW
using Plots
using ProgressMeter

#----

param = Parameters( ϵ  = 1/2,
                    N  = 2^12,
                    L  = 10,
                    T  = 5.0,
                    dt = 0.001)

bump    = Bump(param,1)
solver  = RK4(param)
cheng   = CGBSW(param)
times   = Times(param.dt, param.T)

#----

function create_animation( bump, solver, cheng, times )


    h = cheng.Pi .* fft(bump.h)
    u = cheng.Pi .* fft(bump.u)

    prog = Progress(times.Nt,1)

    hr = real(similar(h))

    anim = @animate for l in range(1,times.Nt-1)

        dt = times.t[l+1]-times.t[l]

        step!(solver, cheng, h, u, dt)

        p = plot(layout=(2,1))

        hr = real(ifft(h))

        plot!(p[1,1], cheng.mesh.x, hr;
	          ylims=(-0.6,1),
        	  title="physical space",
              label=cheng.label)

        plot!(p[2,1], fftshift(cheng.mesh.k),
              log10.(1e-18.+abs.(fftshift(h)));
        	  title="frequency",
          label=cheng.label)

        next!(prog)

    end when mod(l, 200) == 0

    gif(anim, "anim.gif", fps=15); nothing

end

# @time create_animation( bump, solver, cheng, times )

#----
#md # ![](cheng.gif)
