# # Workspace
#
#md # [`notebook`](@__NBVIEWER_ROOT_URL__notebooks/two_problems.ipynb)
#
using WaterWaves1D, Plots, FFTW;

gr();

#---- parameters
param = (
    μ = 0.01,
    ϵ = 0.1,
    N = 2^8, # number of collocation points
    L = π,# size of the mesh (-L,L)
    T = 0.1,# final time of computation
    dt = 0.001,  # timestep
    Ns=1
);

K=floor(param.N/2*2/3)-1
#---- Set initial data
ζ(x,y) = 1/2*cos.(x).*cos.(y)';
ux(x,y) = cos.(y').*sin.(x)+cos.(K*y').*sin.(K*x)/K;
uy(x,y) = -sin.(y').*cos.(x)-sin.(K*y').*cos.(K*x)/K;

init = Init2D(ζ, ux, uy);
# ζ et ux,uy sont des fonctions.
# Init les mets sous une forme utilisable par le programme

#---- Build problems


model=SaintVenant2D(param; dealias = 1 , hamiltonian = false, smooth=0)
modelfast=SaintVenant2D_fast(param; dealias = 1 , hamiltonian = false, smooth=0,large_data=false)
modelfast2=SaintVenant2D_fast(param; dealias = 1 , hamiltonian = false, smooth=0,large_data=true)


solver=RK4(model.mapto(init))


problem = Problem(  model, init, param ; solver=solver )
problem_fast = Problem(  modelfast, init, param ; solver=solver )
problem_fast2 = Problem(  modelfast2, init, param ; solver=solver )

#---- Solve problems

@time solve!(problem_fast)
@time solve!(problem_fast2)
@time solve!(problem)


# Sanity check
η,vx,vy,x,y=problem.model.mapfro(problem.data.U[end])
ηf,vxf,vyf,xf,yf=problem_fast.model.mapfro(problem_fast.data.U[end])
ηf2,vxf2,vyf2,xf2,yf2=problem_fast2.model.mapfro(problem_fast2.data.U[end])

using LinearAlgebra
norm(η-ηf),norm( vx-vxf), norm( vy-vyf)
norm(η-ηf2),norm( vx-vxf2), norm( vy-vyf2)

problem_fast ≈ problem_fast2,problem_fast ≈ problem
nothing
# # Profiling allocations to improve the code
# using BenchmarkTools
# using Profile,PProf

# U=problem_fast.data.U[1]
# step!(solver, modelfast, U, 1)

# Profile.Allocs.clear()
# Profile.Allocs.@profile sample_rate=1 step!(solver, modelfast, U, 1)
# PProf.Allocs.pprof(from_c=false)



#---- Solve other problems



model_hamiltonian=SaintVenant2D_fast(param; dealias = 1, hamiltonian = true)
model_smooth=SaintVenant2D_fast(param; dealias = 1, smooth = true, hamiltonian = false)

problem_hamiltonian = Problem(  model_hamiltonian, init, param ; solver=solver )
problem_smooth = Problem(  model_smooth, init, param ; solver=solver )

@time solve!(problem_hamiltonian)
@time solve!(problem_smooth)

#---- Visualization

η,vx,vy, = solution(problem_fast;T=1)
hη,hvx,hvy, = solution(problem_hamiltonian;T=1)
sη,svx,svy, = solution(problem_smooth;T=1)


mesh=Mesh(param)
kx=mesh.k;ky=mesh.k;x=mesh.x;y=mesh.x;
∂y = 1im * ky' ;
∂x = 1im .* kx;

function ∂fy( f )
    real.(ifft(∂y .* fft(f, 2), 2))
end
function ∂fx( f )
    real.(ifft(∂x .* fft(f, 1), 1))
end


plot(x,y,vx,st=:surface)
plot(x,y,hvx,st=:surface)
plot(x,y,svx,st=:surface)
plot(fftshift(kx),fftshift(ky),log10.(abs.(fftshift(ifft(vx)))),st=:surface)
plot(fftshift(kx),fftshift(ky),log10.(abs.(fftshift(ifft(hvx)))),st=:surface)
plot(fftshift(kx),fftshift(ky),log10.(abs.(fftshift(ifft(svx)))),st=:surface)

plot(x,y,η-0*sη,st=:surface)
plot(x,y,hη,st=:surface)
plot(x,y,sη,st=:surface)
plot(fftshift(kx),fftshift(ky),log10.(abs.(fftshift(ifft(η)))),st=:surface)
plot(fftshift(kx),fftshift(ky),log10.(abs.(fftshift(ifft(hη)))),st=:surface)
plot(fftshift(kx),fftshift(ky),log10.(abs.(fftshift(ifft(sη)))),st=:surface)

plot(x,y,∂fx(vx),st=:surface)

plot(x,y,∂fx(∂fx(vx)),st=:surface)
plot(x,y,∂fx(∂fx(hvx)),st=:surface)
plot(x,y,∂fx(∂fx(svx)),st=:surface)


@gif for time in LinRange(0,1,101)
    η,vx,vy, = solution(problem;T=time)
    plot(x,y,η,st=:surface,zlims=(-1,2));
end
# anim = @animate for time in LinRange(0,param.T,101)
#     η,vx,vy, = sol(problem;T=time)
#     plot(x,y,η,st=:surface);
#     zlims!(-1, 2)
# end
# gif(anim,"SaintVenant2D.gif")

# Energy preservation
function energy(η,vx,vy)
    (η.^2 .+ (1 .+ η) .* (vx.^2 .+ vy.^2))/2
end

function energy_norm(η,vx,vy)
    real.(fft(energy(η,vx,vy)))[1,1]/length(η)
end

function energy_norm(problem,time)
    η,vx,vy, = solution(problem;T=time)
    energy_norm(η,vx,vy)
end

times=LinRange(0,1,11)
plot(times,[energy_norm(problem,t) for t in times])
plot!(times,[energy_norm(problem_hamiltonian,t) for t in times])
plot!(times,[energy_norm(problem_smooth,t) for t in times])