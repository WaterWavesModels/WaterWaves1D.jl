# # Workspace
#
#md # [`notebook`](@__NBVIEWER_ROOT_URL__notebooks/two_problems.ipynb)
#
using WaterWaves1D, Plots, FFTW;

gr();

#---- function to define parameters
function param(; ϵ = 1, N = 2^8, L = 5, T = 1, dt = 1e-3, Ns = 100)
    (
    ϵ = ϵ, # nonlinearity parameter
    N = N, # number of collocation points
    L = L, # size of the mesh (-L,L)
    T = T, # final time of computation
    dt = dt,  # timestep
    Ns = Ns,  # number of stored timesteps
    );
end

#---- define initial data


# function in H^{α+1/2} with one singularity
α = 1.5   # set-up the regularity
ff(x) =  0.5*exp(-abs(x)^α).*exp(-2*abs(x)^2);  


# function in H^{α} with dyadic Fourier contributions
α = 2   # set-up the regularity
ff(x) =  0.5*cos.(x).*exp(-2*abs(x)^2) .+
    0.5*cos.(2^2*x).*exp(-2*abs(x)^2)/2^(2α).+
    0.5*cos.(2^3*x).*exp(-2*abs(x)^2)/2^(3α).+
    0.5*cos.(2^4*x).*exp(-2*abs(x)^2)/2^(4α).+
    0.5*cos.(2^5*x).*exp(-2*abs(x)^2)/2^(5α).+
    0.5*cos.(2^6*x).*exp(-2*abs(x)^2)/2^(6α).+
    0.5*cos.(2^7*x).*exp(-2*abs(x)^2)/2^(7α).+
    0.5*cos.(2^8*x).*exp(-2*abs(x)^2)/2^(8α).+
    0.5*cos.(2^9*x).*exp(-2*abs(x)^2)/2^(9α).+
    0.5*cos.(2^10*x).*exp(-2*abs(x)^2)/2^(10α).+
    0.5*cos.(2^11*x).*exp(-2*abs(x)^2)/2^(11α).+
    0.5*cos.(2^12*x).*exp(-2*abs(x)^2)/2^(12α)


ζ(x) = ff.(x);
u(x) = cos.(x.+1).*ff.(x); #5*(x.+.05).*f.(x.+.5) ;
init = Init(ζ, u);

# # random function in H^s
mesh=Mesh(param(;N=2^18))
f=random(mesh.x,s=2,λ=2,L=2)/2
g=random(mesh.x,s=2,λ=2,L=2)/2
init=Init(mesh,f,g,fast=true)
# init = Random(param(;N=2^18);s=2,λ=2,L=2,a=(1/2,1/2))


#---- plot the initial data
mesh=Mesh(param(;L=5,N=2^15))
plot(mesh.x,[init.η(mesh.x) init.v(mesh.x)],label=["η" "v"])
plot(mesh.k,[abs.(fft(init.η(mesh.x))).+eps() abs.(fft(init.v(mesh.x))).+eps()],label=["η" "v"],yaxis=:log)

#---- function to define problem
function pb( parameters = param() ;  initial_data = init, dealias = 1, label=label)
    Problem( SaintVenant(parameters;  dealias = dealias, label = label) , initial_data, parameters)
end

#---- set problems
problems = Problem[];
Ns = 2 .^(6:15);
for N in Ns
    push!(problems, pb(param(N=N,L=5,dt=1e-4), dealias=1,label="N=2^$(Int(log(2,N)))"))
end

#---- solve problems
for problem in problems
    @time solve!(problem)
end

#---- visualization
plot(problems;T=0,var=[:surface,:fourier])
plot(problems;T=1,var=[:surface,:fourier])

#---- compute norms of errors
using LinearAlgebra

# differentiation function
function diff(u,x) 
    k=Mesh(x).k
    real(ifft(1im*k.*fft(u)))
end

# compute norms
function Compute(;T=1)
    (η,v,x,t)=solution(problems[end],T=T) # reference solution
    p=plot()    # initialize plot
    E=Float64[] # initialize errors
    dE=Float64[] # initialize errors

    for pb in problems[1:end-1]
        label = pb.label
        interpolation_factor =  Int(length(x)/length(pb.data.U[1][:,1]))
        print(interpolation_factor)
        (ηₐ,vₐ)=solution(pb,interpolation=interpolation_factor,T=T) # approximate solution
        plot!(p,x,[abs.(η-ηₐ) .+ eps()  abs.(v-vₐ) .+ eps()],label=["η $label" "v $label"])
        push!(E,norm(norm(η-ηₐ,2)+norm(v-vₐ),2)/norm(norm(η,2)+norm(v),2))
        push!(dE,norm(norm(diff(η,x)-diff(ηₐ,x),2)+norm(diff(v,x)-diff(vₐ,x)),2)/norm(norm(diff(η,x),2)+norm(diff(v,x)),2))
    end
    display(p)
    return (E,dE,p)
end

# Errors at initial times
E₀,dE₀,p₀ = Compute(T=0);
print(log.(E₀[1:end-1]./E₀[2:end])./log(2)) # power of decay (L² norm)
print(log.(dE₀[1:end-1]./dE₀[2:end])./log(2)) # power of decay (̇H¹ norm)

# Errors at final times
E,dE,p = Compute(T=1);
print(log.(E[1:end-1]./E[2:end])./log(2)) # power of decay (L² norm)
print(log.(dE[1:end-1]./dE[2:end])./log(2)) # power of decay (̇H¹ norm)

# Plots
scatter(Ns[1:end-1],E[1:end],xaxis=:log,yaxis=:log, label = "T=1")
scatter!(Ns[1:end-1],E₀[1:end],xaxis=:log,yaxis=:log, label = "T=0")
plot!(Ns[1:end-1],100*Ns[1:end-1].^(-2),label="N^-2")

scatter!(Ns[1:end-1],dE[1:end],xaxis=:log,yaxis=:log, label = "T=1")
scatter!(Ns[1:end-1],dE₀[1:end],xaxis=:log,yaxis=:log, label = "T=0")
plot!(Ns[1:end-1],10*Ns[1:end-1].^(-1),label="N^-1")

title!("L² and ̇H¹ norms of the error at initial and final time")
xlabel!("N")
nothing

########## 
# This concerns the Burgers equations
##########

# include("../src/models/Burgers.jl")


# param = (
#     ϵ = 1,
#     N = 2^14, # number of collocation points
#     L = 10,# size of the mesh (-L,L)
#     T = .1,# final time of computation
#     dt = 0.0001,  # timestep
#     Ns = 1000, # number of stored timesteps
# );


# #---- initial data
# fh(x) =  0.5*exp(-abs(x*5)^3);


# ζ(x) = fh.(x);
# u(x) = 5*(x.+.05).*fh.(x.+.5) ;
# init = Init(ζ, u);

# pb1= Problem( Burgers(param;  dealias = 1) , init, param, solver=RK4(param,1))
# pb2= Problem( Burgers(param;  dealias = 10) , init, param)
# pb3= Problem( Burgers(param;  dealias = 50) , init, param)
# pb4= Problem( Burgers(param;  dealias = 100) , init, param)


# problems = [pb1 pb2 pb3 pb4] 

# #---- computation
# for problem in problems
#     @time solve!(problem)
# end

# #---- visualization
# plot(problems;T=0,var=[:surface,:fourier])
# plot(problems;T=0.1,var=[:surface])
# plot(problems;T=0.1,var=[:fourier])
# plot!(problems;T=0,var=[:fourier])

# xlims!(-10,10)