using DeepWaterModels
using FFTW
using Plots

function fig_problem!( p, problem::Problem )

    s = 0
    (hhat,uhat) = problem.data[end]
    (hr,ur)     = (real(ifft((problem.model.Gamma.^s).*hhat)),
		   real(ifft(uhat)))
    
    plot!(p[1,1], problem.model.mesh.x, hr; 
		  title="physical space",
	          label=problem.model.label)
    plot!(p[2,1], fftshift(problem.model.mesh.k),
                  log10.(1e-18.+abs.(fftshift(hhat))); 
		  title="frequency",
    	          label=problem.model.label)  
    
end

param = Parameters( ϵ  = 1/2, 
                    N  = 2^12,
                    L  = 10,
                    T  = 5,
                    dt = 0.001)

bump     = Bump( param )
solver   = RK4( param.N )

times    = Times(param.dt, param.T)

cheng    = Cheng( bump.mesh, param.ϵ)
problem1 = Problem( cheng, bump, param, solver )

matsuno  = Matsuno(bump.mesh, param.ϵ)
problem2 = Problem(matsuno, bump, param, solver )

solve!( problem1, times )
solve!( problem2, times )

p = plot(layout=(2,1))

fig_problem!( p, problem1 )
fig_problem!( p, problem2 )

display(p)
