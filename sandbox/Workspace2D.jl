# # Workspace
#
#md # [`notebook`](@__NBVIEWER_ROOT_URL__notebooks/two_problems.ipynb)
#
using WaterWaves1D, Plots, FFTW;

gr();

#---- parameters
param = (
    μ = 0.01,
    ϵ = 1,
    N = 2^9, # number of collocation points
    L = π,# size of the mesh (-L,L)
    T = 0.5,# final time of computation
    dt = 0.001,  # timestep
    Ns=10
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
model=SaintVenant2D(param; dealias = 1 , hamiltonian = false)
model_hamiltonian=SaintVenant2D(param; dealias = 1, hamiltonian = true)
model_smooth=SaintVenant2D(param; dealias = 1/50, hamiltonian = false)

solver=RK4(model.mapto(init))

problem = Problem(  model, init, param ; solver=solver )
problem_hamiltonian = Problem(  model_hamiltonian, init, param ; solver=solver )
problem_smooth = Problem(  model_smooth, init, param ; solver=solver )


#---- Solve problems
solve!(problem)
solve!(problem_hamiltonian)
solve!(problem_smooth)


#---- Visualization


function sol(pb;T=nothing)
    if isnothing(T) T = pb.times.ts[end] end
	T=min(max(T,0),pb.times.ts[end])
	index = findfirst(pb.times.ts.>=T)
	t=pb.times.ts[index]
    return (pb.model.mapfro)(pb.data.U[index])
end


η,vx,vy, = sol(problem;T=1)
hη,hvx,hvy, = sol(problem_hamiltonian;T=1)
sη,svx,svy, = sol(problem_smooth;T=1)


mesh=Mesh(param)
kx=mesh.k;ky=mesh.k;x=mesh.x;y=mesh.x;
∂y = 1im * ky' .* one.( x);
∂x = one.(y') .* 1im .* kx;

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
    η,vx,vy, = sol(problem;T=time)
    plot(x,y,η,st=:surface,zlims=(-1,2));
end
# anim = @animate for time in LinRange(0,param.T,101)
#     η,vx,vy, = sol(problem;T=time)
#     plot(x,y,η,st=:surface);
#     zlims!(-1, 2)
# end
# gif(anim,"SaintVenant2D.gif")

function energy(η,vx,vy)
    (η.^2 .+ (1 .+ η) .* (vx.^2 .+ vy.^2))/2
end

function energy_norm(η,vx,vy)
    real.(fft(energy(η,vx,vy)))[1,1]/length(η)
end

function energy_norm(problem,time)
    η,vx,vy, = sol(problem;T=time)
    energy_norm(η,vx,vy)
end

times=LinRange(0,1,11)
plot(times,[energy_norm(problem,t) for t in times])
plot!(times,[energy_norm(problem_hamiltonian,t) for t in times])
plot!(times,[energy_norm(problem_smooth,t) for t in times])