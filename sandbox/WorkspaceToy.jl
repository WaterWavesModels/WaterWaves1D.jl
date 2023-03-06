# # Workspace
#
#md # [`notebook`](@__NBVIEWER_ROOT_URL__notebooks/two_problems.ipynb)
#
using WaterWaves1D, Plots, FFTW;

gr();

#---- parameters
param = (
    δ = 0.1,
    a = 100,
    N = 2^10, # number of collocation points
    L = π,# size of the mesh (-L,L)
    T = 1,# final time of computation
    dt = 0.0001,  # timestep
);


#---- initial data
ζ(x) = 1 .+ sin.(x)/2;
u(x) = -sin.(x);
init = Init(ζ, u);
# ζ et v sont des fonctions.
# Init les mets sous une forme utilisable par le programme


#---- models to compare
#models = AbstractModel[]


pb= Problem( Toy(param; dealias = 1) , init, param)


#problems = Problem[]
#for model in models
#    push!(problems, Problem(model, init, param;solver=s))
#end
problems = [pb] 
#problems = [pWW pGN pFGa pFGb pFGc] 

#---- computation
for problem in problems
    @time solve!(problem)
end

#---- visualization
plt=plot(pb;T=1,var=[:surface,:velocity])
plt=plot!(pb;T=0.0,var=[:surface,:velocity])



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
    u = ifft(pb.data.U[index][:,2])
    h = pb.data.U[index][:,1]
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