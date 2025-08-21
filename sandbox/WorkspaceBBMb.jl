# # Workspace
#
using WaterWaves1D, Plots, FFTW;

gr();

include("./models/BBM.jl")
include("./models/HBBM.jl")
include("./models/HBBMb.jl")
include("./models/KG.jl")

function solution_full(pb;T=nothing)
    if isnothing(T) T = pb.times.ts[end] end
	T=min(max(T,0),pb.times.ts[end])
	index = findfirst(pb.times.ts.>=T)
	t=pb.times.ts[index]
	
    return pb.model.mapfrofull(pb.data.U[index])
end

function solution_BBM(pb;δ,T=nothing)
    if isnothing(T) T = pb.times.ts[end] end
	T=min(max(T,0),pb.times.ts[end])
	index = findfirst(pb.times.ts.>=T)
	t=pb.times.ts[index]

		
    _,u,x=pb.model.mapfro(pb.data.U[index])
    ∂ₓ	=  1im * Mesh(x).k
	F = ∂ₓ./(1 .+ abs.(δ * ∂ₓ ).^2)

    w=real.(ifft(∂ₓ.*fft(u)) )
    v=real.(ifft(∂ₓ.*F.*fft(u .+ 1/2*u.^2)))
	
    return u,v,w,x
end



#---- parameters
param = (
    δ = 0.01,
    ϵ = 0.01,
    N = 2^8, # number of collocation points
    L = 5,# size of the mesh (-L,L)
    T = 0.5,# final time of computation
    dt = 1e-4,  # timestep
    Ns = 1000, # number of stored timesteps
);



#---- initial data
fu(x) =  0.5*exp(-(x)^2);

function pu(x)
    h = 0. ;
    for i in -100:100
        h += fu(x + 5*i ) 
    end
    return h
end

u₀(x) = fu.(x);
b(x)=sin.(2π*x)/4

x=Mesh(param).x
plot(x,[u₀(x) b(x)])

#---- models to compare
#models = AbstractModel[]

init = Init(b,u₀)

pBBM= Problem( BBM(param; dealias = 1) , init, param)
pBBM0= Problem( HBBMb(param; dealias = 1, prep=1) , init, param)
pBBM1= Problem( HBBMb(param; dealias = 1, prep=1,stable=0) , init, param)
pBBM2= Problem( HBBMb(param; dealias = 1, prep=1,stable=1) , init, param)


problems = [pBBM pBBM0 pBBM1 pBBM2] 

#---- computation
for problem in problems
    @time solve!(problem)
end

#---- visualization

uₒ,v₀,w₀,B,x₀=solution_full(pBBM0;T=0.5);
u₁,v₁,w₁,B,x₁=solution_full(pBBM1;T=0.5);
u₂,v₂,w₂,B,x₂=solution_full(pBBM2;T=0.5);

_,uB,X=solution(pBBM;T=0.5);


plot(x,uB-uₒ)
plot!(x,uB-u₁)
plot!(x,uB-u₂)


∂=1im*Mesh(x).k
function diff(f;n=1)
    real(ifft((∂.^n).*fft(f)))
end


plot(x,diff(uB-uₒ;n=4))
plot!(x,diff(uB-u₁;n=4))
plot!(x,diff(uB-u₂;n=4))

plot(x,[v₀ v₁ v₂])
plot(x,[w₀ w₁ w₂])


k=Mesh(x).k;δ=param.δ;
F = 1 ./(1 .+ (δ*k).^2);
uᵢ = real.(ifft(F.*(fft(u)-δ^2*1im*k.*fft(w))))
plot(x,[uB-u,uB-uᵢ] )
plot(x,[uB-u₂,uB-uᵢ] )

plot(x,uᵢ)


param = (
    δ = 0.01,
    ϵ = 0.0001,
    N = 2^11, # number of collocation points
    L = 100,# size of the mesh (-L,L)
    T = 0.1,# final time of computation
    dt = 1e-5,  # timestep
    Ns = 1000, # number of stored timesteps
);


pKG = Problem( KG(param; dealias = 1, c₀ = 1) , Init(x->-x.*u₀(x),zero), param)
solve!(pKG)

plot(pKG,T=0.01,var=[:surface,:velocity])

v,w,x=solution(pKG,T=0.1)
maximum(abs.(v))
maximum(abs.(w))



uB,vB,wB,xB=solution_BBM(pBBM;δ=param.δ,T=0.1)
u,v,w,x=solution_full(pBBM2;T=0.1)
plot(x,[u-uB (v-vB)*param.ϵ*param.δ (w-wB)*param.δ],label=["u" "ϵδv" "δw"])

anim = @animate for time in LinRange(0,param.T,501)
    plot((pBBM,pBBM0); T = time,var=[:velocity])
    ylims!(-3e-4, 3e-4)
end
gif(anim,"BBM0.gif")


anim = @animate for time in LinRange(0,param.T,501)
    plot((pBBM,pBBM1); T = time,var=[:velocity])
    ylims!(-1e-5, 1e-5)
end
gif(anim,"BBM1.gif")

anim = @animate for time in LinRange(0,param.T,501)
    plot((pBBM,pBBM2); T = time,var=[:velocity])
    ylims!(-1e-5, 1e-5)
end
gif(anim,"BBM2.gif")

"""
    Norm(problem,kwargs)


Returns the `W^{k,p}` norm of the m-th time-derivative of the function `u` at time `t=T``
## Optional keyword arguments
- parameter `T` (by default last computed time)
- parameters `(m,k,p)` (by default `(0,0,Inf)`)
- `dot` (default = `false`): returns the homogeneous Sobolev norm.
"""
function norm(x,u;k=0,p=Inf,dot=false)
    ∂ₓ=1im*Mesh(x).k
    
    if dot == false
        Js=0:k
    else
        Js=k
    end 
    N=0
    for j in Js
        if p == Inf
            N += maximum(abs.(ifft( ( ∂ₓ).^j .* fft(u) )))
        else
            N += (sum( abs.(ifft( ( ∂ₓ).^j .* fft(u) )).^p )*Mesh(x).dx/2π).^(1/p)
        end
    end
    return N
end


function norm(pb::Problem;T=Inf::AbstractFloat,k=0,p=Inf,dot=false)
    _,u,x=solution(pb;T=T);    
    return norm(x,u;k=k,p=p,dot=dot)
end

function norm(pb1::Problem,pb2::Problem;T=Inf::AbstractFloat,k=0,p=Inf,dot=false)
    _,u1,_=solution(pb1;T=T);    
    _,u2,x=solution(pb2;T=T);    

    return norm(x,u1-u2;k=k,p=p,dot=dot)
end

norm(pBBM,pBBM2;T=1,p=2)

function normfast(pb::Problem;T=Inf::AbstractFloat,k=0,p=Inf,dot=false)
    u,v,w,x=solution_full(pb;T=T);    
    return norm(x,v;k=k,p=p,dot=dot),norm(x,w;k=k,p=p,dot=dot)
end

function Diff( ; δ, ϵ, prep=1, dt = 1e-4, T=0.5, L = 10, N = 2^9, p=2)
    param = (
        δ = δ,
        ϵ = ϵ,
        N = N, # number of collocation points
        L = L,# size of the mesh (-L,L)
        T = T,# final time of computation
        dt = dt,  # timestep
        Ns = 1, # number of stored timesteps
    );


    #---- initial data
    u₀(x) = 0.5*exp.(-(x).^2);

    x=Mesh(param).x
    
    #---- models to compare
    #models = AbstractModel[]

    init = Init(zero,u₀)

    pBBM= Problem( BBM(param; dealias = 1) , init, param)
    phBBM= Problem( HBBM(param; dealias = 1, prep=prep) , init, param)
    
    #---- computation
    solve!(pBBM)
    solve!(phBBM)
    
    return norm(pBBM,phBBM;p=p)

end

N=[]
δ=1
prep = 1
eps= 10 .^(-(1:0.5:4))
for ϵ in eps
    push!(N,Diff(δ=δ,ϵ=ϵ,prep=prep,dt=1e-5,L=100,N=2^11,p=Inf))
end

scatter(eps,N,xscale=:log10,yscale=:log10,label="difference (prep = $prep)")
if prep == 1
    plot!(eps,δ^2*eps,label="\$δ^2 ϵ\$")
    plot!(eps,eps.^(2),label="\$ϵ^{2}\$")
elseif prep == 2
    plot!(eps,δ^2*eps.^(2),label="\$δ^2ϵ^{2}\$")
end
    plot!(eps,eps.^(2),label="\$ϵ^{2}\$")

N=[]
del= 10 .^(-(0:0.25:2))
for δ in del
    push!(N,Diff(δ=δ,ϵ=0.01,prep=2))
end

scatter(del,N,xscale=:log10,yscale=:log10,label="difference (prep = $prep)")
plot!(del,0.0001*del.^(2),label="\$δ^2ϵ\$")



function Normfast( ; δ, ϵ, prep=1, dt = 1e-4, T=0.5, L = 10, N = 2^9, p=2)
    param = (
        δ = δ,
        ϵ = ϵ,
        N = N, # number of collocation points
        L = L,# size of the mesh (-L,L)
        T = T,# final time of computation
        dt = dt,  # timestep
        Ns = 1, # number of stored timesteps
    );


    #---- initial data
    u₀(x) = 0.5*exp.(-(x).^2);

    x=Mesh(param).x
    
    #---- models to compare
    #models = AbstractModel[]

    init = Init(zero,u₀)
    pBBM= Problem( BBM(param; dealias = 1) , init, param)
    phBBM= Problem( HBBM(param; dealias = 1, prep=prep) , init, param)
    
    #---- computation
    solve!(pBBM)
    solve!(phBBM)
    
    _,u,_=solution(pBBM;T=T);    
    _,v,w,x=solution_full(phBBM;T=T);    
    k=Mesh(x).k
    du = real.(ifft(1im*k.*fft(u)))
    return norm(x,w.-du;k=0,p=p)

end



N=[]
δ=0.1
prep = 1
eps= 10 .^(-(1:0.5:4))
for ϵ in eps
    push!(N,Normfast(δ=δ,ϵ=ϵ,prep=prep,dt=1e-5,T=0.1))
end


scatter!(eps,[n[1] for n in N],xscale=:log10,yscale=:log10,label="difference (prep = $prep)")
if prep == 1
    plot!(eps,δ*eps,label="\$δ ϵ\$")
    plot!(eps,eps.^(3/2),label="\$ϵ^{3/2}\$")
elseif prep == 2
    plot!(eps,δ^2*eps.^(2),label="\$δ^2ϵ^{2}\$")
end