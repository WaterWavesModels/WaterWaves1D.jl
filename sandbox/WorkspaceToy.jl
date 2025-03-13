# # Workspace
#
#md # [`notebook`](@__NBVIEWER_ROOT_URL__notebooks/two_problems.ipynb)
#
using WaterWaves1D, Plots, FFTW;

gr();

#---- parameters
param = (
    δ = 0.025,
    a = 1000,
    N = 2^10, # number of collocation points
    L = 10,# size of the mesh (-L,L)
    T = 1,# final time of computation
    dt = 0.0001,  # timestep
    Ns = 1000, # number of stored timesteps
);

param2 = (
    δ = 0.025,
    a = 10000,
    N = 2^10, # number of collocation points
    L = 10,# size of the mesh (-L,L)
    T = 1,# final time of computation
    dt = 0.00001,  # timestep
    Ns = 1000, # number of stored timesteps
);


#---- initial data
fh(x) =  0.5*exp(-(x)^2);

function Fh(x)
    h = 1. ;
    for i in -100:100
        h += fh(x + 5*i ) 
    end
    return h
end

ζ(x) = 1 .+ fh.(x.+5);
u(x) = (x.+5).*exp.(-(x.+5).^2) / (1 .+ 0* param.a);
init = Init(ζ, u);
# ζ et v sont des fonctions.
# Init les mets sous une forme utilisable par le programme


#---- models to compare
#models = AbstractModel[]


pb1= Problem( Toy(param; str = true , dealias = 0) , init, param)
pb2= Problem( Toy(param2; str = true , dealias = 0) , init, param2)


problems = [pb1 pb2] 

#---- computation
for problem in problems
    @time solve!(problem)
end

#---- visualization
plot(problems;T=1,var=[:velocity])

anim = @animate for time in LinRange(0,param.T,501)
    plot(problems; T = time,var=[:velocity])
    ylims!(-0.5, 0.5)
end
gif(anim,"toy.gif")



function energy(pb,param;T=Inf::AbstractFloat)
    mesh=Mesh(param)
    T=min(max(T,0),pb.times.ts[end])
	index = findfirst(pb.times.ts.>=T)
    u = ifft(pb.data.U[index][2])
    h = pb.data.U[index][1]
    h,u = pb.model.mapfrofull(pb.data.U[index])
    
    return sum(abs.(u).^2 .*h)
end
E1=[];E2=[]
for t in Times(param).ts
    push!(E1,energy(pb1,param;T=t))
    push!(E2,energy(pb2,param;T=t))
end
plot(Times(param).ts,[E1 E2])

# Attention au scaling a revoir !
# Faire la derivee en temps
"""
    Norm(problem,kwargs)


Returns the `W^{k,p}` norm of the m-th time-derivative of the function `u` at time `t=T``
## Optional keyword arguments
- parameter `T` (by default last computed time)
- parameters `(m,k,p)` (by default `(0,0,Inf)`)
- `dot` (default = `false`): returns the homogeneous Sobolev norm.
- `δ` (default = `1`): returns the scaled Sobolev norm with prefactors associated with `u(δ⋅)`.
"""
function norm(pb,param;T=Inf::AbstractFloat,k=0,p=Inf,m=0,δ=1,dot=false)

    mesh=Mesh(param)
    ∂ₓ=1im*mesh.k
    Jδ = param.δ * ∂ₓ .+ 1im

    T=min(max(T,0),pb.times.ts[end])
	index = findfirst(pb.times.ts.>=T)
    u = ifft(pb.data.U[index][2])
    h = pb.data.U[index][1]
    h,u = pb.model.mapfrofull(pb.data.U[index])
    # h,u = solution(pb,T=T) #pb, you have only the real part

    for n in 1:m
        u = -param.a * ifft(Jδ  .* fft(u) ) ./ h 
    end
    if dot == false
        Js=0:k
    else
        Js=k
    end 
    N=0
    for j in Js
        if p == Inf
            N += maximum(abs.(ifft( (δ * ∂ₓ).^j .* fft(u) )))
        else
            N += (sum( abs.(ifft( (δ * ∂ₓ).^j .* fft(u) )).^p )*mesh.dx/2π).^(1/p)
        end
    end
    return N
end


N=[]
for t in Times(param).ts
    push!(N,norm(pb,param;T=t,m=1,k=0,p=Inf,dot=true))
end
plot(Times(param).ts,N)


# savefig("err.pdf")