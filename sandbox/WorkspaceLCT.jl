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
    N = 2^10, # number of collocation points
    L = π, # size of the mesh (-L,L)
    T = 1, # final time of computation
    dt = 0.001,  # timestep
);


#---- initial data
ζ0(x) = zero(x);
u0(x) = -sin.(x);
init = Init(ζ0, u0);
# ζ et v sont des fonctions.
# Init les mets sous une forme utilisable par le programme


pWW= Problem( WaterWaves(param; dealias = 1) , init, param)
pGN= Problem( WhithamGreenNaghdi(param; SGN=true, dealias = 1) , init, param)
pLCT=Problem[]
as = [20 40 60 80 100 120 140 160 180 200]
for a in as
    para=merge(param,(a=a,))
    push!(pLCT , Problem(relaxedGreenNaghdi(para;iterate=false,FG=false, id = 0, dealias = 1),init,para) )
end
solve!([pWW pGN])
solve!(pLCT)

plt=plot([pWW pGN pLCT],var=[:surface,:velocity])

@gif for time in LinRange(0,Times(param).tfin,201)
    plt=plot([pWW,pGN,pLCT],var=[:surface,:velocity],T=time,legend=:topright,xlims=(0,3),ylims=(-1,2));
    #plot!(plt,[pFG3b,pGN],var=[:differences],T=time,legend=:bottomleft,xlims=(0,10),ylims=(-1e-5,1e-5));
end




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

    T=min(max(T,0),pb.times.ts[end])
	index = findfirst(pb.times.ts.>=T)
    
    U = copy(pb.data.U[index])
    Um = copy(pb.data.U[index-1])
    (η,v,u,q,w,x) = pb.model.mapfrofull(U)
    (ηm,vm,um,pm,wm,xm) = pb.model.mapfrofull(Um)

    pb.model.f!(U)
    pb.model.f!(Um)

    (dtη,~,dtu,dtq,dtw,~) = pb.model.mapfrofull(U)
    (dtηm,~,dtum,dtqm,dtwm,~) = pb.model.mapfrofull(Um)

    (d2tη,d2tu,d2tq,d2tw)=(dtη-dtηm,dtu-dtum,dtq-dtqm,dtw-dtwm)./pb.times.dt


    function |(u,δ,j,p)
        if p == Inf
            maximum(abs.(ifft( (δ * ∂ₓ).^j .* fft(u) )))
        else
            (sum( abs.(ifft( (δ * ∂ₓ).^j .* fft(u) )).^p )*mesh.dx/2π).^(1/p)
        end
    end


    if dot == false
        Js=0:k
    else
        Js=k
    end 
    N=0
    if m == 0
        for j in Js
            N += |(η,δ,j,p) + |(u,δ,j,p) 
        end
    elseif m == 1 
        for j in Js[1:end-1]
            N += |(dtη,δ,j,p) + |(dtu,δ,j,p) 
        end
    elseif m == 2 
        for j in Js[1:end-2]
            N += |(d2tη,δ,j,p) + |(d2tu,δ,j,p) 
        end
    end
    return copy(N)
end

N=[]
for pb in pLCT
    push!(N,norm(pb,param,T=0.25))
end

#---- models to compare
#models = AbstractModel[]
function Compare(scenario;values=nothing)
    if scenario == 1
        @info("different values")

    end
end

# pFGb = Problem(relaxedGreenNaghdi(p1;FG=false, id = 1, dealias = 1),init,p1,solver=RK4_naive(),label="a=20")  # dealias = 1 pour dealiasing
# pFGc = Problem(relaxedGreenNaghdi(p1;FG=false, id = 2, dealias = 1),init,p1,solver=RK4_naive(),label="a=20")  # dealias = 1 pour dealiasing

# pFG1a = Problem(relaxedGreenNaghdi(p1;iterate=false,id=0, dealias = 1),init,p1,solver=RK4_naive(),label="a=20 (id0)")  # dealias = 1 pour dealiasing
# pFG1b = Problem(relaxedGreenNaghdi(p1;iterate=true,id=1, dealias = 1),init,p1,solver=RK4_naive(),label="a=20 (id1)")  # dealias = 1 pour dealiasing
# pFG1c = Problem(relaxedGreenNaghdi(p1;id=2, dealias = 1),init,p1,solver=RK4_naive(),label="a=20 (id2)")  # dealias = 1 pour dealiasing

# p2=merge(param,(a=100,))
# pFG2a = Problem(relaxedGreenNaghdi(p2;id=0, dealias = 1),init,p1,solver=RK4_naive(),label="a=100 (id0)")  # dealias = 1 pour dealiasing
# pFG2b = Problem(relaxedGreenNaghdi(p2;id=1, dealias = 1),init,p1,solver=RK4_naive(),label="a=100 (id1)")  # dealias = 1 pour dealiasing
# pFG2c = Problem(relaxedGreenNaghdi(p2;id=2, dealias = 1),init,p1,solver=RK4_naive(),label="a=100 (id2)")  # dealias = 1 pour dealiasing

# p3=merge(param,(a=200,))
# pFG3a = Problem(relaxedGreenNaghdi(p3;id=0, dealias = 1),init,p1,solver=RK4_naive(),label="a=200 (id0)")  # dealias = 1 pour dealiasing
# pFG3b = Problem(relaxedGreenNaghdi(p3;id=1, dealias = 1),init,p1,solver=RK4_naive(),label="a=200 (id1)")  # dealias = 1 pour dealiasing
# pFG3c = Problem(relaxedGreenNaghdi(p3;id=2, dealias = 1),init,p1,solver=RK4_naive(),label="a=200 (id2)")  # dealias = 1 pour dealiasing


#problems = Problem[]
#for model in models
#    push!(problems, Problem(model, init, param;solver=s))
#end
# problems = [pWW pGN pFG1a pFG2a pFG3a pFG1b pFG2b pFG3b pFG1c pFG2c pFG3c] 
#problems = [pWW pGN pFGa pFGb pFGc] 

#---- computation
# for problem in problems
#     @time solve!(problem)
# end

#---- visualization
plt=plot(problems[1];T=1,var=[:surface])
plt=plot(problems[2];T=1,var=[:surface,:velocity,:fourier])
plt=plot!(problems[3];T=1,var=[:surface,:velocity,:fourier])
plt=plot!(problems[4];T=1,var=[:surface,:velocity,:fourier])
plt=plot!(problems[5];T=1,var=[:surface,:velocity,:fourier])

plt=plot(problems[2];T=1,var=[:surface,:velocity,:fourier])
plt=plot!(problems[6];T=1,var=[:surface,:velocity,:fourier])
plt=plot!(problems[7];T=1,var=[:surface,:velocity,:fourier])
plt=plot!(problems[8];T=1,var=[:surface,:velocity,:fourier])

plt=plot(problems[2];T=1,var=[:surface,:velocity,:fourier])
plt=plot!(problems[9];T=1,var=[:surface,:velocity,:fourier])
plt=plot!(problems[10];T=1,var=[:surface,:velocity,:fourier])
plt=plot!(problems[11];T=1,var=[:surface,:velocity,:fourier])

#plot(problems,var=[:differences],T=2)

plot(problems[[2;3]],var=[:difference,:difference_velocity],T=1,legend=:bottomleft)
plot(problems[[2;6]],var=[:difference,:difference_velocity],T=1,legend=:bottomleft)
plot(problems[[2;9]],var=[:difference,:difference_velocity],T=1,legend=:bottomleft)
plot!(problems[[2;4]],var=[:difference,:difference_velocity],T=1,legend=:bottomleft)
plot!(problems[[2;7]],var=[:difference,:difference_velocity],T=1,legend=:bottomleft)
plot!(problems[[2;10]],var=[:difference,:difference_velocity],T=1,legend=:bottomleft)
plot!(problems[[2;5]],var=[:difference,:difference_velocity],T=1,legend=:bottomleft)
plot(problems[[2;8]],var=[:difference,:difference_velocity],T=1,legend=:bottomleft)
plot!(problems[[2;11]],var=[:difference,:difference_velocity],T=1,legend=:bottomleft)


plot!(problems[[1;2]],var=[:difference,:difference_velocity],T=1,legend=:bottomleft)


plot(problems[[3;4]],var=[:differences])

plot(problems[[2;3;6]],var=[:differences],T=2,legend=:bottomleft)

plot(problems[[2;3]],var=[:differences],T=2)
plot!(problems[[2;4]],var=[:differences],T=2)

plot(problems[[2;4]],var=[:differences],T=2)
plot!(problems[[2;5]],var=[:differences],T=2)



@gif for time in LinRange(0,2,101)
    plt=plot([pFG3a,pGN],var=[:differences],T=time,legend=:bottomleft,xlims=(0,10),ylims=(-1e-5,1e-5));
    plot!(plt,[pFG3b,pGN],var=[:differences],T=time,legend=:bottomleft,xlims=(0,10),ylims=(-1e-5,1e-5));

end

# Attention au scaling a revoir !
# Faire la derivee en temps



function norm(p;T=Inf::AbstractFloat,a,μ,s=0)
    T=min(max(T,0),p.times.ts[end])
	index = findfirst(p.times.ts.>=T)
	η,v,u,p,w,x = (p.model.mapfrofull)(p.data.U[index])
	mesh=Mesh(x)
    ∂ₓ=1im*mesh.k
    |(x)=maximum(abs.(x))
    
    if s==0
        return |(η),|(u),|(p),|(w)
    else
        return |(ifft(∂ₓ.^s.*fft(η))),|(ifft(∂ₓ.^s.*fft(u))),|(ifft(∂ₓ.^s.*fft(p))),|(ifft(∂ₓ.^s.*fft(w)))
    end

end

function norms(p;T=nothing,a,μ,s=0,rel=false)
    if T == nothing T = p.times.ts end
    N=zeros(length(T),4)
    j=0
    if rel == true N0 = sum(norm(p;T=T[1],a=a,μ=μ,s=s)) else N0=1 end

    for t in T
        j+=1
        n=norm(p;T=t,a=a,μ=μ,s=s)
        N[j,:].=n./N0
    end
    return N
end

function norm(p1,p2;T=Inf::AbstractFloat,a,μ)
    T=min(max(T,0),p1.times.ts[end],p2.times.ts[end])
	i1 = findfirst(p1.times.ts.>=T)
	η1,v1,u1 = (p1.model.mapfrofull)(p1.data.U[i1])
	i2 = findfirst(p2.times.ts.>=T)
	η2,v2,u2 = (p2.model.mapfrofull)(p2.data.U[i2])
	mesh=Mesh(x)
    ∂ₓ=1im*mesh.k
    |(x)=maximum(abs.(x))
    
    return |(η1-η2),|(u1-u2)
    
end

function norms(p1,p2;T=nothing,a,μ,rel=false)
    if T == nothing T = p1.times.ts end
    N=zeros(length(T),2)
    j=0
    if rel == true N0 = sum(norm(p1,p2;T=T[1],a=a,μ=μ)) else N0=1 end

    for t in T
        j+=1
        n=norm(p1,p2;T=t,a=a,μ=μ)
        N[j,:].=n./N0
    end
    return N
end


plot(norms(pFG3b;a=100,μ=0.01,s=2),label=["η" "u" "p" "w"])

plot(norms(pGN,pFG1c;a=100,μ=0.01),label=["Δη" "Δu"])


errid0=[sum(norm(pGN,pFG1a;a=100,μ=0.01)),
    sum(norm(pGN,pFG2a;a=100,μ=0.01)),
    sum(norm(pGN,pFG3a;a=100,μ=0.01))]
errid1=[sum(norm(pGN,pFG1b;a=100,μ=0.01)),
    sum(norm(pGN,pFG2b;a=100,μ=0.01)),
    sum(norm(pGN,pFG3b;a=100,μ=0.01))]
errid2=[sum(norm(pGN,pFG1c;a=100,μ=0.01)),
    sum(norm(pGN,pFG2c;a=100,μ=0.01)),
    sum(norm(pGN,pFG3c;a=100,μ=0.01))]


scatter([20,100,200],errid0,label="no preparation")
scatter!([20,100,200],errid1,label="w well-prepared")
scatter!([20,100,200],errid2,label="w and p well-prepared")
plot!(xscale=:log10,yscale=:log10,legend=:bottomleft,
        xlabel="a",ylabel="difference (L∞)",
        title="difference between GN and LCT")


scatter!([(200,sum(norm(pGN,pFG3c;a=100,μ=0.01)))])
scatter!([(20,sum(norm(pGN,pFG1c;a=20,μ=0.01)[1]))])
norm(pGN,pFG2c;a=100,μ=0.01)

savefig("err.pdf")