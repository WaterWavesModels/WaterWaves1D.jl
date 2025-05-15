# # Workspace
#
#md # [`notebook`](@__NBVIEWER_ROOT_URL__notebooks/two_problems.ipynb)
#
using WaterWaves1D,FFTW,Plots,LinearAlgebra,ProgressMeter;

gr();

#---- parameters
param = (
    ϵ = 1,
    N = 2^9, # number of collocation points
    L = π,# size of the mesh (-L,L)
    T = 0.1,# final time of computation
    dt=1e-3,  # timestep
    Ns=1
);

K=floor(param.N/3)

using LinearAlgebra

mesh=Mesh(param)
kx=mesh.k;ky=mesh.k;x=mesh.x;y=mesh.x;
∂y = 1im * ky' ;
∂x = 1im * kx  ;

function ∂fy( f )
    real.(ifft(∂y .* fft(f, 2), 2))
end
function ∂fx( f )
    real.(ifft(∂x .* fft(f, 1), 1))
end

#-----
#Experiment 0. Compare 2D with 1D  (I coded this before noticing you did the same)

#---- Initialize problems
ζ(x,y) = 0.5 * (cos.(x).+cos.(K*x)/K).*ones(length(y))';
ux(x,y) = 0.5 * sin.(x).*ones(length(y))';
uy(x,y) = zero(x).*ones(length(y))';

init = Init(x->ζ(x,0),x-> ux(x,0));
init2D = Init2D(ζ, ux, uy);

# I also checked with ζ,uy depending only on y, ux=0
#ζ(x,y) = 0.5 * (cos.(y).+cos.(K*y)/K)'.*ones(length(x));
#uy(x,y) = 0.5 * sin.(y)'.*ones(length(x));
#ux(x,y) = zero(x).*ones(length(y))';
#init = Init(x->ζ(0,x)'[:],x-> uy(0,x)'[:]);
#init2D = Init2D(ζ, ux, uy);


# I checked dealias ∈ {0,1} , hamiltonian ∈ {false,true}, smooth ∈ {0,1}
model = SaintVenant(param; dealias = 1, smooth = 0)
model2D=SaintVenant2D_fast(param; dealias = 1 , hamiltonian = false, smooth=0)

solver=RK4(model.mapto(init))
solver2D=RK4(model2D.mapto(init2D))

problem = Problem(  model, init, param ; solver=solver )
problem2D = Problem(  model2D, init2D, param ; solver=solver2D )

#---- Solve problems
solve!(problem)
solve!(problem2D)

η,v= solution(problem);
η2D,vx,vy,x,y = solution(problem2D);

#---- Perform checks
η2D[:,1] ≈ η
sum(η2D,dims=2) ≈ param.N* η
vx[:,1] ≈ v
sum(vx,dims=2) ≈ param.N* v
vy==zero(vy)
sum(abs.(η2D[:,1]- η).+abs.(vx[:,1]- v)+abs.(vy[:,1]))

plot(x,[η2D[:,1], η])
xlims!(-0.1,0.1)
ylims!((maximum(η)-0.01,maximum(η)+0.01))

#-----
#Experiment 0.5 Compare hamiltonian with non-hamiltonian

#---- Initialize problems

# This initial data is irrotational. Change coefficients or signs to violate irrotationality
ζ(x,y) = 0.5 .*cos.(x).*cos.(y)';
ux(x,y) = 0.5 .* cos.(y').*sin.(x)+cos.(K*y').*sin.(K*x)/K;
uy(x,y) = 0.5 .* sin.(y').*cos.(x)+sin.(K*y').*cos.(K*x)/K;

init = Init2D(ζ, ux, uy);

# I checked dealias ∈ {0,1} , smooth ∈ {0,1}
# When dealias = 0 (also when dealias = 3/4) there is a significant difference ; maybe some spurious instabilities from numerical errors ?
model=SaintVenant2D_fast(param; dealias = 1/4 , hamiltonian = false, smooth=1)
model_hamiltonian=SaintVenant2D_fast(param; dealias = 1/4 , hamiltonian = true, smooth=1)

solver=RK4(model.mapto(init))

problem = Problem(  model, init, param ; solver=solver )
problem_hamiltonian = Problem(  model_hamiltonian, init, param ; solver=solver )

#---- Solve problems
@time solve!(problem)
@time solve!(problem_hamiltonian)

#---- Perform checks
η,vx,vy,x,y = solution(problem;T=1)
hη,hvx,hvy,hx,hy = solution(problem_hamiltonian;T=1)

hη ≈ η
hvx ≈ vx
hvy ≈ vy
sum(abs.(hη - η).+abs.(hvx - vx)+abs.(hvy - vy))

plot(x,y,hη ,st=:surface)
plot!(x,y,η ,st=:surface)
plot(x,y,η- hη ,st=:surface)

# In order to explore the issue when dealias < 1
plot(x,y,fftshift(abs.(fft(η))) ,st=:surface,scale=:log10)
plot(x,y,fftshift(abs.(fft(hη))) ,st=:surface,scale=:log10)
plot(x,y,fftshift(abs.(fft(hη- η))) ,st=:surface,scale=:log10)

plot(x[1:200],y[1:200],fftshift(abs.(fft(η)))[1:200,1:200] ,st=:surface,scale=:log10)
plot(x[1:200],y[1:200],fftshift(abs.(fft(hη)))[1:200,1:200] ,st=:surface,scale=:log10)
plot(x[1:200],y[1:200],fftshift(abs.(fft(hη- η)))[1:200,1:200] ,st=:surface,scale=:log10)

maximum(fftshift(abs.(fft(η)))[1:200,1:200])
maximum(fftshift(abs.(fft(hη)))[1:200,1:200])
maximum(fftshift(abs.(fft(η-hη)))[1:200,1:200])
# Both models suffer from dealiasing instabilities


#-----
#Experiment 1. Zero-deth. Sharp cut-off. Initial data in H^2.  
#Clear instabilities if negative depth, but not so clear at zero depth.

#---- Initialize problems
ζ(x,y) = -1.05 .*cos.(x).*cos.(y)';
ux(x,y) = 0.5 .* cos.(y').*sin.(x)+cos.(K*y').*sin.(K*x)/K^2;
uy(x,y) = -0.5 .* sin.(y').*cos.(x)+ 2 .* sin.(K*y').* cos.(K*x)/K^2;

init = Init2D(ζ, ux, uy);


model=SaintVenant2D_fast(param; dealias = 1 , hamiltonian = false, smooth=0)
model_hamiltonian=SaintVenant2D_fast(param; dealias = 1 , hamiltonian = true, smooth=0)

#smooth models for comparison
model_smooth=SaintVenant2D_fast(param; dealias = 1 , hamiltonian = false, smooth=1)
model_hamiltonian_smooth=SaintVenant2D_fast(param; dealias = 1 , hamiltonian = true, smooth=1)

solver=RK4(model.mapto(init))

problem = Problem(  model, init, param ; solver=solver )
problem_hamiltonian = Problem(  model_hamiltonian, init, param ; solver=solver )
problem_smooth = Problem(  model_smooth, init, param ; solver=solver )
problem_hamiltonian_smooth = Problem(  model_hamiltonian_smooth, init, param ; solver=solver )

#---- Solve problems
@time solve!(problem)
@time solve!(problem_hamiltonian)
@time solve!(problem_smooth)
@time solve!(problem_hamiltonian_smooth)

#---- Visualization

η,vx,vy,x,y = solution(problem;T=1)
hη,hvx,hvy,hx,hy = solution(problem_hamiltonian;T=1)
ηs,vxs,vys,xs,ys = solution(problem_smooth;T=1)
hηs,hvxs,hvys,hxs,hys = solution(problem_hamiltonian_smooth;T=1)

plot(x,y,∂fx(∂fx(vx)),st=:surface,label="non-hamilt.")
plot(x,y,∂fx(∂fx(vxs)),st=:surface,label="non-hamilt. smooth")
plot(x,y,∂fx(∂fx(hvx)),st=:surface, label="hamilt.")
plot(x,y,∂fx(∂fx(hvxs)),st=:surface, label="hamilt. smooth")

nothing



#Experiment 2. Violating 1 + \eta + abs(v)^2 > 0. Sharp cut-off. Initial data in H^2. 
#Decreasing or increasing N greatly affects only hamiltonian system with sharp cut-off.
#This reflects the ill-posedness of the hamiltonian system. Presumably for N very large even with smooth cut-off we should see instabilities.
#---- Initialize problems
ζ(x,y) = -0.5*cos.(x).*cos.(y');
ux(x,y) = 2 .*cos.(y').*sin.(x)+cos.(K*y').*sin.(K*x)/K^2;
uy(x,y) = - 2 .* sin.(y').*cos.(x);

init = Init2D(ζ, ux, uy);

model=SaintVenant2D_fast(param; dealias = 1 , hamiltonian = false, smooth=0)
model_hamiltonian=SaintVenant2D_fast(param; dealias = 1 , hamiltonian = true, smooth=0)

#smooth models for comparison
model_smooth=SaintVenant2D_fast(param; dealias = 1 , hamiltonian = false, smooth=1)
model_hamiltonian_smooth=SaintVenant2D_fast(param; dealias = 1 , hamiltonian = true, smooth=1)

solver=RK4(model.mapto(init))

problem = Problem(  model, init, param ; solver=solver )
problem_hamiltonian = Problem(  model_hamiltonian, init, param ; solver=solver )
problem_smooth = Problem(  model_smooth, init, param ; solver=solver )
problem_hamiltonian_smooth = Problem(  model_hamiltonian_smooth, init, param ; solver=solver )

#---- Solve problems
@time solve!(problem)
@time solve!(problem_hamiltonian)
@time solve!(problem_smooth)
@time solve!(problem_hamiltonian_smooth)

#---- Visualization

η,vx,vy,x,y = solution(problem;T=1)
hη,hvx,hvy,hx,hy = solution(problem_hamiltonian;T=1)
ηs,vxs,vys,xs,ys = solution(problem_smooth;T=1)
hηs,hvxs,hvys,hxs,hys = solution(problem_hamiltonian_smooth;T=1)

plot(x,y,∂fx(∂fx(vx)),st=:surface,label="non-hamilt.")
plot(x,y,∂fx(∂fx(vxs)),st=:surface,label="non-hamilt. smooth")
plot(x,y,∂fx(∂fx(hvx)),st=:surface, label="hamilt.")
plot(x,y,∂fx(∂fx(hvxs)),st=:surface, label="hamilt. smooth")

nothing


#Experiment 4: Zero-depth. y-independent problem. Comparison of 2D and 1D solver. 
#---- Initialize problems
f(x) = 1;
ζ(x,y) = -1*cos.(x).*f.(y)';
ux(x,y) = 2 .* f.(y').*sin.(x)+ f.(y').*sin.(K*x)/K^2;
uy(x,y) = 0 .* f.(y').*f.(x);

init1D = Init(x->-1*cos.(x), x->2*sin.(x).+sin.(K*x)/K^2);
init2D = Init2D(ζ, ux, uy);


model1D = SaintVenant_fast(param; dealias = true, smooth = false, label = "1D") 
model2D=SaintVenant2D_fast(param; dealias = 1 , hamiltonian = true, smooth=0)

solver1D=RK4(model1D.mapto(init1D))
solver2D=RK4(model2D.mapto(init2D))

problem1D = Problem(model1D, init1D, param;solver=solver1D)
problem2D = Problem(model2D, init2D, param ; solver=solver2D)

#---- Solve problems
solve!(problem1D)
@time solve!(problem2D)


#---- Visualization

η2D,vx2D,vy2D,x2D,y2D = solution(problem2D;T=0.1)
η1D,v1D,x1D = solution(problem1D;T=0.1)

plot(x2D,y2D,vx2D,st=:surface)
plot(x1D,v1D)

∂=1im*Mesh(x).k
diff(η;n=1)=real.(ifft(∂.^n.*fft(η)))

plot(x2D,y2D,∂fx(∂fx(vx2D)),st=:surface)
plot(x1D,diff(v1D;n=2))

plot(x1D,[diff(v1D;n=2) ∂fx(∂fx(vx2D))[:,1]],label=["1D" "2D"])
nothing



# I tried to make convergence plots, but this part(below) is not working. 

#Convergence for experiment 1 with optional minimal depth. 
function IntegrateSV(;ϵ=1,L=π,N=2^9,T=0.1,dt =1e-3,dealias=1,smooth=false,Ns=nothing,hamiltonian = false)
	if isnothing(Ns)
		params = ( ϵ  = ϵ,
				N  = N, L  = L,
	            T  = T, dt = dt )
	else
		params = ( ϵ  = ϵ,
				N  = N, L  = L,
	            T  = T, dt = dt,
				Ns = Ns )
	end


    #K=floor(params.N/2*2/3)-1 # If data depends on K, then the reference solution needs to be rebuilt for each value of K
    K = 16;
    ζ(x,y) = 0.5 .*cos.(x).*cos.(y)' .+ 0*cos.(K*y').*sin.(K*x)/K^2; 
    ux(x,y) = 1 .* cos.(y').*sin.(x) .+ 0*cos.(K*y').*sin.(K*x)/K^2; 
    uy(x,y) = -sin.(y').*cos.(x);

    init = Init2D(ζ, ux, uy);

    model=SaintVenant2D_fast(params; dealias = dealias , hamiltonian = hamiltonian, smooth=smooth)
    solver=RK4(model.mapto(init))
    problem = Problem(  model, init, params ; solver=solver )
    @time solve!(problem)

	return problem
end


function Convergence(smooth=false,hamiltonian = false)
	problems=Problem[]
	E0=Float64[]
	E1=Float64[]

    Nref=2^11
    reference_problem=IntegrateSV(L=π,N=Nref,T=0.1,dt = 1e-4,smooth=smooth,hamiltonian=hamiltonian)
    ηref,vxref, vyref,xref, yref=solution(reference_problem; T= 0.05)
    Uref=reference_problem.data.U[end]/Nref^2
    kxref = fftshift(Mesh(xref).k)
    Uref1 = fftshift(Uref[1])
    Uref2 = fftshift(Uref[2])
    Uref3 = fftshift(Uref[3])

    for n=6:10
        p = IntegrateSV(L=π,N=2^n,T=0.1,dt = 1e-4,Ns=1,smooth=smooth,hamiltonian=hamiltonian)
        push!(problems,p)
        η,vx, vy,x, y=solution(p; T= 0.05)
        #U=[p.data.U[end][1] ; p.data.U[end][2] ; p.data.U[end][3]]
        U = p.data.U[end]/2^(2n)
        U1 = fftshift(U[1])
        U2 = fftshift(U[2])
        U3 = fftshift(U[3])
        Nmodes=length(x)÷2
        Ucomp1 = Uref1[Nref÷2-Nmodes+1:Nref÷2+Nmodes, Nref÷2-Nmodes+1:Nref÷2+Nmodes];
        Ucomp2 = Uref2[Nref÷2-Nmodes+1:Nref÷2+Nmodes, Nref÷2-Nmodes+1:Nref÷2+Nmodes];
        Ucomp3 = Uref3[Nref÷2-Nmodes+1:Nref÷2+Nmodes, Nref÷2-Nmodes+1:Nref÷2+Nmodes];

        kx = fftshift(Mesh(x).k)
        println(kx)
        println(size(Ucomp1))
        println(size(U1))

        pl=plot(x, y,η,st=:surface)
        display(pl)

        # pl=plot(kx, kx,log10.(abs.(U1)),st=:surface)
        # display(pl)
        # plot!(kx,kx,log10.(abs.(Ucomp1)),st=:surface)
        # display(pl)

        #pl2=plot(fftshift(kx),fftshift(kx),log10.(abs.(fftshift(Ucomp1))),st=:surface)
        #pl2=plot(kx, kx,log10.(abs.(Ucomp1)),st=:surface)
        #pl3=plot(kx,kx,log10.(abs.(U1-Ucomp1)),st=:surface)
        #display(pl2)
        #display(pl3)

        push!(E0,norm(norm(U1-Ucomp1,2)+norm(U2-Ucomp2,2)+ norm(U3 - Ucomp3,3),2)/norm(norm(Ucomp1,2)+norm(Ucomp2, 2)+norm(Ucomp3,2),2))

        # k=[sqrt.(1 .+ Mesh(x).k.^2);sqrt.(1 .+ Mesh(x).k.^2)]
        # push!(E1,norm(k.*U-k.*Ucomp,2)/norm(k.*Ucomp,2))
    end
	return reference_problem,E0 #,E1
end


#--- Figures
plot_font = "Computer Modern"
default(fontfamily=plot_font)

#--- Experiment 1 : heap of water, sharp low-pass filter
# Solve the problems, save the errors  and the reference solution
reference_problem,E0=Convergence()

# Plot convergence rates
Fig1b=plot(;xlabel="\$N\$",axis=:log)
Ns=2 .^(6:10);
#scatter!(Fig1b,Ns,E1,label="\$E_1\$",color=1)
scatter!(Fig1b,Ns,E0,label="\$E_0\$",color=2)
#plot!(Fig1b,Ns,Ns.^(-1),label="",color=1)
plot!(Fig1b,Ns,100 .*Ns.^(-2),label="",color=2)

savefig(Fig1b,"Fig1b.pdf");savefig(Fig1b,"Fig1b.svg");

