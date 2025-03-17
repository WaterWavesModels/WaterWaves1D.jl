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
    T = 1,# final time of computation
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


# Profiling
using BenchmarkTools
using Profile,PProf

U=problem_fast.data.U[1]
step!(solver, modelfast, U, 1)

Profile.Allocs.clear()
Profile.Allocs.@profile sample_rate=1 step!(solver, modelfast, U, 1)
PProf.Allocs.pprof(from_c=false)


model_hamiltonian=SaintVenant2D(param; dealias = 1, hamiltonian = true)
model_smooth=SaintVenant2D(param; dealias = 1/50, hamiltonian = false)



problem_hamiltonian = Problem(  model_hamiltonian, init, param ; solver=solver )
problem_smooth = Problem(  model_smooth, init, param ; solver=solver )




f  = rand(ComplexF64,(2^12,2^12));
fᵗ = similar(f);
a = similar(f);
b = similar(f);
c = similar(f);

FFTW.set_num_threads(4)
Px = plan_fft(f,  1, )#flags=FFTW.PATIENT)    
Py = plan_fft(fᵗ, 1, )#flags=FFTW.PATIENT)
P2 = plan_fft(f,  2, )#flags=FFTW.PATIENT)    

iPx = plan_ifft(f,  1)#, flags=FFTW.PATIENT)    
iPy = plan_ifft(fᵗ, 1)# flags=FFTW.PATIENT)
iP2 = plan_ifft(f,  2)#, flags=FFTW.PATIENT)    


using LinearAlgebra
function my_ifft!(fft,f,fᵗ,gᵗ)
    mul!(f, iPx, fft )
    transpose!(fᵗ,f)
    mul!(gᵗ, iPy, fᵗ )
    transpose!(fft,gᵗ)
end

function my_ifft2!(fft,f)
    mul!(f, iPx, fft )
    mul!(fft, iP2, f )
end
function my_ifft3!(fft,f,fᵗ,gᵗ)
    ldiv!(f, Px, fft )
    transpose!(fᵗ,f)
    ldiv!(gᵗ, Py, fᵗ )
    transpose!(fft,gᵗ)
end

function my_ifft4!(fft,f)
    ldiv!(f, Px, fft )
    ldiv!(fft, P2, f )
end

function my_ifft5!(fft)
    ifft!(fft)
end


my_ifft!(f,a,b,c);
my_ifft2!(f,fᵗ);
my_ifft3!(f,a,b,c);
my_ifft4!(f,fᵗ);
my_ifft5!(f);
ifft!(f);


@info("0")
display(@benchmark ifft!(f))

@info("1")
display(@benchmark my_ifft!(f,a,b,c))
@info("2")
display(@benchmark my_ifft2!(f,a))
@info("3")
display(@benchmark my_ifft3!(f,a,b,c))
@info("4")
display(@benchmark my_ifft4!(f,a))
@info("5")
display(@benchmark my_ifft5!(f,a))



function with_transpose(;Nsteps=2^6,Nmesh=2^6)
f  = rand(ComplexF64,(Nmesh,Nmesh))
f̂  = similar(f)
fᵗ = rand(ComplexF64,(Nmesh,Nmesh))
f̂ᵗ = similar(fᵗ)
    
FFTW.set_num_threads(2)
Px = plan_fft(f,  1)# flags=FFTW.PATIENT)    
Py = plan_fft(fᵗ, 1)# flags=FFTW.PATIENT)
    
    
    for n=1:Nsteps     
        ldiv!(f, Px, f̂ )
        transpose!(f̂ᵗ,f)
        ldiv!(fᵗ, Py, f̂ᵗ)
        transpose!(f,fᵗ)        
    end

    real(f)

end
function without_transpose(;Nsteps=2^6,Nmesh=2^6)
    f  = rand(ComplexF64,(Nmesh,Nmesh))
    f̂  = similar(f)
    fᵗ = rand(ComplexF64,(Nmesh,Nmesh))
        
    FFTW.set_num_threads(4)
    Px = plan_fft(f,  1)# flags=FFTW.PATIENT)    
    Py = plan_fft(fᵗ, 2)# flags=FFTW.PATIENT)
        
        
        for n=1:Nsteps     
            ldiv!(fᵗ, Px, f̂ )
            ldiv!(f, Py, fᵗ)
        end
    
        real(f)
    
    end
    

with_transpose(Nmesh=2^8,Nsteps=2^10)
without_transpose(Nmesh=2^8,Nsteps=2^10)

with_fft_transposed_bench = @benchmark with_transpose(Nmesh=2^8,Nsteps=2^10)
without_fft_transposed_bench = @benchmark without_transpose(Nmesh=2^8,Nsteps=2^10)


#display(@benchmark mul!(f, iPx, a))
#display(@benchmark mul!(f, iPy, a))
#display(@benchmark mul!(f, iP2, a))


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