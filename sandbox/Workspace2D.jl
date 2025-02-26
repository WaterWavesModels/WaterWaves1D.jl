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
    N = 2^8, # number of collocation points
    L = π,# size of the mesh (-L,L)
    T = 1,# final time of computation
    dt = 0.01,  # timestep
);


#---- initial data
ζ(x,y) = exp.(-x.*x).*exp.(-y.*y)';
ux(x,y) = cos.(y').*sin.(x);
uy(x,y) = sin.(y').*cos.(x);

init = Init2D(ζ, ux, uy);
# ζ et ux,uy sont des fonctions.
# Init les mets sous une forme utilisable par le programme

# mesh=Mesh(param)
# kx=mesh.k;ky=mesh.k;x=mesh.x;y=mesh.x;
# ∂y = 1im * ky' .* one.( x)
# ∂x = one.(y') .* 1im .* kx

# function ∂fy( f )
#     ifft(∂y .* fft(f, 2), 2)
# end
# function ∂fx( f )
#     ifft(∂x .* fft(f, 1), 1)
# end

# f = sin.(2*x) .* cos.(3*y') 
# g = exp.(-x.*x.*one.(y')-one.(x).*y'.*y')

# plotly()

# plot(x,y,f,st = [:surface, :contourf])
# plot(x,y,real.(∂fx(f)),st = [:surface, :contourf])

# real.(∂fx(f))-2*cos.(2*x).*cos.(3*y)'

# plot(x,y,g, st = :surface)
# plot(x,y,real.(∂fx(g)), st = :surface)

# plot(real.(∂fy(g))+20*y'.*exp.(-10*x.*x).*exp.(-10*y.*y)', st = :surface)
# g=exp.(-10*x.*x).*exp.(-10*y.*y)'

#---- models to compare
#models = AbstractModel[]

model=SaintVenant2D(param; dealias = 1 , hamiltonian = false)
model_hamiltonian=SaintVenant2D(param; dealias = 1, hamiltonian = true)
model_smooth=SaintVenant2D(param; dealias = 1/50, hamiltonian = false)

solver=RK4(Saint_Venant_model.mapto(init))
problem = Problem(  model, init, param ; solver=solver )
problem_hamiltonian = Problem(  model_hamiltonian, init, param ; solver=solver )
problem_smooth = Problem(  model_smooth, init, param ; solver=solver )



solve!(problem)
solve!(problem_hamiltonian)
solve!(problem_smooth)


#---- visualization


function sol(pb;T=nothing)
    if isnothing(T) T = pb.times.ts[end] end
	T=min(max(T,0),pb.times.ts[end])
	index = findfirst(pb.times.ts.>=T)
	t=pb.times.ts[index]
    return (pb.model.mapfro)(pb.data.U[index])
end


η,vx,vy, = sol(problem)
hη,hvx,hvy, = sol(problem_hamiltonian)
sη,svx,svy, = sol(problem_smooth)


mesh=Mesh(param)
kx=mesh.k;ky=mesh.k;x=mesh.x;y=mesh.x;
plot(x,y,η,st=:surface)
plot(x,y,hη,st=:surface)
plot(fftshift(kx),fftshift(ky),log10.(abs.(fftshift(ifft(η)))),st=:surface)
plot(fftshift(kx),fftshift(ky),log10.(abs.(fftshift(ifft(hη)))),st=:surface)

plot(x,y,η-sη,st=:surface)
plot(fftshift(kx),fftshift(ky),log10.(abs.(fftshift(ifft(η)))),st=:surface)
plot(fftshift(kx),fftshift(ky),log10.(abs.(fftshift(ifft(sη)))),st=:surface)


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