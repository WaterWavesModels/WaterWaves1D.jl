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
pyplot()

#----

param = Parameters( ϵ  = 1/2, 
                    N  = 2^12,
                    L  = 10,
                    T  = 10.0,
                    dt = 0.001)

bump    = Bump(param)
solver  = RK4(param)
cheng   = Cheng(param)
times   = Times(param.dt, param.T)

#----

function create_animation( bump, solver, cheng, times )

    N  = cheng.mesh.N
    h  = zeros(ComplexF64, N)
    u  = zeros(ComplexF64, N)
    
    bump(h, u)
    
    h .= cheng.Pi .* fft(h)
    u .= cheng.Pi .* fft(u)
                   
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
    
    end when mod(l, 100) == 0
    
    gif(anim, "cheng.gif", fps=15); nothing 

end

@time create_animation( bump, solver, cheng, times )

#----
#md # ![](cheng.gif)
