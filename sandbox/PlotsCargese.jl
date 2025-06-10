# # Workspace
#
using WaterWaves1D, FFTW, LinearAlgebra, Plots;
# using Plots; gr();


param = (
    ϵ = 1,
    N = 2^10, # number of collocation points
    L = π, # size of the mesh (-L,L)
    T = 1, # final time of computation
    dt = 0.0001,  # timestep
    );

#---- initial data
ζ0(x) = zero(x);
u0(x) = -sin.(x);
init = Init(ζ0, u0);
# ζ0 et u0 sont des fonctions.
# Init les mets sous une forme utilisable par le programme


#pWW= Problem( WaterWaves(param; dealias = 1) , init, param)
#pGN= Problem( WhithamGreenNaghdi(param; SGN=true, dealias = 1) , init, param)
#solve!([pWW pGN])


as = [20 100 500]

pLCT=Problem[]
for a in as
    para=merge( param,( μ = 0.01, a=a ) )
    push!(pLCT , Problem(relaxedGreenNaghdi(para; id = 1, dealias = 1),init,para) )
end

solve!(pLCT)
plot(pLCT[1],T=0,var=[:surface,:velocity],label="")
plot(pLCT,label=["ε=0.05" "ε=0.01" "ε=0.002"])


η1,=solution(pLCT[1])
η2,=solution(pLCT[2])
η3,=solution(pLCT[3])

mesh=Mesh(param)
∂=1im*mesh.k
d4η1=real.(ifft((∂.^4).*fft(η1)))
d4η2=real.(ifft((∂.^4).*fft(η2)))
d4η3=real.(ifft((∂.^4).*fft(η3)))

plot(mesh.x,[d4η1 d4η2 d4η3],label=["ε=0.05" "ε=0.01" "ε=0.002"],
    title="4th derivative of the surface deformation",
    titlefontsize = 10,
    xlabel="x",
    ylabel="∂⁴η ")
#xlims!(0,2)

d8η1=real.(ifft((∂.^8).*fft(η1)))
d8η2=real.(ifft((∂.^8).*fft(η2)))
d8η3=real.(ifft((∂.^8).*fft(η3)))

plot(mesh.x,[d8η1 d8η2 d8η3],label=["ε=0.05" "ε=0.01" "ε=0.002"],
    title="8th derivative of the surface deformation",
    titlefontsize = 10,
    xlabel="x",
    ylabel="∂⁸η "
    )
#xlims!(0,2)
    
"""
    Wkp(x,u,k,p)


Compute the `L^p` norm of the j-th space-derivative of the function 
whose values at collocation points `x` are `u`.
"""
function Wkp(x,u,k,p)

    if k == 0
        ∂ₓ = one.(x)
        dx = (maximum(x)-minimum(x))/(length(x)-1)
    else
        mesh=Mesh(x)
        ∂ₓ=1im*mesh.k
        dx = mesh.dx
    end

    norm(ifft( ∂ₓ.^k .* fft(u) ) , p )*(dx/2π)^(1/p)

end

"""
    Norm(problem,kwargs)


Return the `W^{k-m,p}` norm of the m-th time-derivative of the function `u` at time `t=T``
## Optional keyword arguments
- parameter `T` (by default last computed time)
- parameters `(m,k,p)` (by default `(0,0,Inf)`)
- `dot` (default = `false`): returns the homogeneous Sobolev norm.
- `δ` (default = `1`): returns the scaled Sobolev norm with prefactors associated with `u(δ⋅)`.
"""
function Norm(pb;T=Inf::AbstractFloat,k=0,p=Inf,m=0,δ=1,dot=false)
    T=min(max(T,0),pb.times.ts[end])
	index = findfirst(pb.times.ts.>=T)
    if dot == false
        Js=0:k
    else
        Js=k
    end 
    N=0
    if m == 0
        U = copy(pb.data.U[index])
        (η,v,u,q,w,x) = pb.model.mapfrofull(U)
        for j in Js
            N += δ^j*( Wkp(x,η,j,p) + Wkp(x,u,j,p) )
        end
    elseif m == 1 
        U = copy(pb.data.U[index])
        pb.model.f!(U)
        (dtη,~,dtu,dtq,dtw,x) = pb.model.mapfrofull(U)
        for j in Js[1:end-1]
            N += δ^j*( Wkp(x,dtη,j,p) + Wkp(x,dtu,j,p) )
        end
    elseif m == 2 
        U = copy(pb.data.U[index])
        pb.model.f!(U)
        (dtη,~,dtu,dtq,dtw,x) = pb.model.mapfrofull(U)
        Um = copy(pb.data.U[index-1])
        pb.model.f!(Um)
        (dtηm,~,dtum,dtqm,dtwm,~) = pb.model.mapfrofull(Um)
        (d2tη,d2tu,d2tq,d2tw)=(dtη-dtηm,dtu-dtum,dtq-dtqm,dtw-dtwm)./pb.times.dt
        for j in Js[1:end-2]
            N += δ^j*( Wkp(x,d2tη,j,p) + Wkp(x,d2tu,j,p) )
        end
    end
    return N
end

# """
#     Norm(pb1,pb2,kwargs)


# Return the `W^{k,p}` norm of the difference at time `t=T``
# ## Optional keyword arguments
# - parameter `T` (by default last computed time)
# - parameters `(k,p)` (by default `(0,Inf)`)
# - `dot` (default = `false`): returns the homogeneous Sobolev norm.
# - `δ` (default = `1`): returns the scaled Sobolev norm with prefactors associated with `u(δ⋅)`.
# """
# function Norm(pb1,pb2;T=Inf::AbstractFloat,k=0,p=Inf,m=0,δ=1,dot=false)
#     η1,v1,x1,t1 = solution(pb1,T=T)
#     η2, = solution(pb2,T=t1,x=x1)
#     η0,v0,x, = solution(pb2,T=0)

#     if dot == false
#         Js=0:k
#     else
#         Js=k
#     end 
#     N=0

#     for j in Js
#         N += δ^j*|(x,η1-η2,j,p) 
#     end
#     return N
# end


"""
    Compute(scenario;WP=1,p=2,k=3)

Compute the W^{k-j,p} norm of the j-th time derivative of the LCT solution, for j in {0,1,2}.

If 'scenario == 1', then the norms are computed for several values of 'a'.

If 'scenario == 2', then the norms are computed for several values of 'δ'.

Optional keyword arguments are 'p' (default = 2), 'k' (default = 3),  
'prep' the preparation of the initial data (default = 1), and
'str' the structre of the system (default = faulse).

"""
function Compute(scenario;prep=1,str=false,p=2,k=3)

    #---- parameters
    param = (
        ϵ = 1,
        N = 2^10, # number of collocation points
        L = π, # size of the mesh (-L,L)
        T = 0.1, # final time of computation
        dt = 0.0001,  # timestep
        Ns=1 #times stored
        );

    #---- initial data
    ζ0(x) = zero(x);
    u0(x) = -sin.(x);
    init = Init(ζ0, u0);
    # ζ0 et u0 sont des fonctions.
    # Init les mets sous une forme utilisable par le programme

    pLCT=Problem[]

    #pWW= Problem( WaterWaves(param; dealias = 1) , init, param)
    #pGN= Problem( WhithamGreenNaghdi(param; SGN=true, dealias = 1) , init, param)
    #solve!([pWW pGN])
    if scenario == 1
        
        as = [20 40 60 80 100 120 140 160 180 200]
        @info("different values of a: $as.\nδ = 0.1")

        for a in as
            para=merge( param,( δ = 0.1, a=a ) )
            push!(pLCT , Problem(threescale(para;str=str, prep = prep, dealias = 1),init,para) )
        end

        solve!(pLCT)

    end

    if scenario == 2

        deltas = 1 ./ sqrt.([10 25 50 75 100 125 150 250 500 1000])
        @info("different values of δ: $deltas.\na=100")

        for delta in deltas
            para=merge( param,( μ = delta^2 , a = 100) )
            push!(pLCT , Problem(relaxedGreenNaghdi(para;FG = str, id = prep, dealias = 1),init,para) )
        end

        solve!(pLCT)

    end

    N = zeros( length(pLCT) , 10, 3 )
    for j in eachindex(pLCT)
        for m in [0 1 2]
            for k in 0:9
                N[j,k+1,m+1]=Norm(pLCT[j],p=p,k=k,m=m)
            end
        end
    end
    return as,N
end

using Plots
as,N=Compute(2;prep=1,str=false)
scatter(as,N[9,:,1]',xaxis=:log10,yaxis=:log10)
nothing
# scatter(mus,NWP3,xaxis=:log10,yaxis=:log10,label=["WW vs SV" "WW vs LCT" "WW vs GN" "LCT vs GN"],legend=:bottomright,title="WP3, a = 100",xlabel="mu",ylabel="difference (elevation, in \$L^\\infty\$)")
# plot!(mus,mus.^2,label="μ²",color=:3)


# """
#     Compare(scenario;WP=1,p=Inf,k=0)

# Compute the W^{k-j,p} norm of the j-th time derivative of the LCT solution, for j in {0,1,2}.

# If 'scenario == 1', then the norms are computed for several values of 'a'.

# If 'scenario == 2', then the norms are computed for several values of 'δ'.

# Optional keyword arguments are 'p' (default = Inf), 'k' (default = 0), and 
# 'WP' the preparation of the initial data (default = 1).
# """
# function Compare(scenario;WP=1,p=Inf,k=0)

#     #---- parameters
#     param = (
#         ϵ = 1,
#         N = 2^10, # number of collocation points
#         L = π, # size of the mesh (-L,L)
#         T = 0.5, # final time of computation
#         dt = 0.001,  # timestep
#         );

#     #---- initial data
#     ζ0(x) = zero(x);
#     u0(x) = -sin.(x);
#     init = Init(ζ0, u0);
#     # ζ0 et u0 sont des fonctions.
#     # Init les mets sous une forme utilisable par le programme

    
#     if scenario == 1
        
#         as = [20 40 60 80 100 120 140 160 180 200]
#         @info("different values of a: $as.\nδ = 0.1")
        
#         para = merge( param,( μ = 0.01,  ) )
#         pWW = Problem( WaterWaves(para; dealias = 1) , init, para) 
#         pGN = Problem( WhithamGreenNaghdi(para; SGN=true, dealias = 1) , init, para) 
#         pSV = Problem( SaintVenant(para; dealias = 1) , init, para) 

#         solve!(pWW);solve!(pGN);solve!(pSV);
        
#         pLCT=Problem[]
#         for a in as
#             para=merge( param,( dt = 0.0001, μ = 0.01, a=a ) )
#             push!(pLCT , Problem(relaxedGreenNaghdi(para;FG=false, id = WP-1, dealias = 1),init,para) )
#         end
#         solve!(pLCT)

#         N = zeros( length(pLCT), 4 )
#         for j in 1:length(pLCT)
#             N[j,1]=Norm(pWW,pSV,p=p,k=k)
#             N[j,2]=Norm(pWW,pLCT[j],p=p,k=k)
#             N[j,3]=Norm(pWW,pGN,p=p,k=k)
#             N[j,4]=Norm(pGN,pLCT[j],p=p,k=k)
#         end
        

#     end

#     if scenario == 2

#         mus = 1 ./ [10 25 50 75 100 125 150 250 500 1000]
#         @info("different values of δ²: $mus.\na=100")


#         pWW=Problem[]
#         pSV=Problem[]
#         pGN=Problem[]
#         pLCT=Problem[]

#         for mu in mus
#             para=merge( param,( dt = 0.0001, μ = mu , a = 100) )
#             push!(pLCT , Problem(relaxedGreenNaghdi(para;FG=false, id = WP-1, dealias = 1),init,para) )
#             para=merge( param,( dt = 0.001, μ = mu , a = 100) )
#             push!(pWW ,  Problem( WaterWaves(para; dealias = 1) , init, para) )
#             push!(pGN ,  Problem( WhithamGreenNaghdi(para; SGN=true, dealias = 1) , init, para) )
#             push!(pSV ,  Problem( SaintVenant(para; dealias = 1) , init, para) )

#         end

#         solve!(pLCT)
#         solve!(pGN)
#         solve!(pWW)
#         solve!(pSV)


#         N = zeros( length(pLCT), 4 )
#         for j in 1:length(pLCT)
#             N[j,1]=Norm(pWW[j],pSV[j],p=p,k=k)
#             N[j,2]=Norm(pWW[j],pLCT[j],p=p,k=k)
#             N[j,3]=Norm(pWW[j],pGN[j],p=p,k=k)
#             N[j,4]=Norm(pGN[j],pLCT[j],p=p,k=k)
#         end

#     end

#     return N
# end




# scatter(mus,NWP3,xaxis=:log10,yaxis=:log10,label=["WW vs SV" "WW vs LCT" "WW vs GN" "LCT vs GN"],legend=:bottomright,title="WP3, a = 100",xlabel="mu",ylabel="difference (elevation, in \$L^\\infty\$)")
# plot!(mus,mus.^2,label="μ²",color=:3)



# plt=plot([pWW pGN pLCT],var=[:surface,:velocity])

# @gif for time in LinRange(0,Times(param).tfin,201)
#     plt=plot([pWW,pGN,pLCT],var=[:surface,:velocity],T=time,legend=:topright,xlims=(0,3),ylims=(-1,2));
#     #plot!(plt,[pFG3b,pGN],var=[:differences],T=time,legend=:bottomleft,xlims=(0,10),ylims=(-1e-5,1e-5));
# end



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
# plt=plot(problems[1];T=1,var=[:surface])
# plt=plot(problems[2];T=1,var=[:surface,:velocity,:fourier])
# plt=plot!(problems[3];T=1,var=[:surface,:velocity,:fourier])
# plt=plot!(problems[4];T=1,var=[:surface,:velocity,:fourier])
# plt=plot!(problems[5];T=1,var=[:surface,:velocity,:fourier])

# plt=plot(problems[2];T=1,var=[:surface,:velocity,:fourier])
# plt=plot!(problems[6];T=1,var=[:surface,:velocity,:fourier])
# plt=plot!(problems[7];T=1,var=[:surface,:velocity,:fourier])
# plt=plot!(problems[8];T=1,var=[:surface,:velocity,:fourier])

# plt=plot(problems[2];T=1,var=[:surface,:velocity,:fourier])
# plt=plot!(problems[9];T=1,var=[:surface,:velocity,:fourier])
# plt=plot!(problems[10];T=1,var=[:surface,:velocity,:fourier])
# plt=plot!(problems[11];T=1,var=[:surface,:velocity,:fourier])

# #plot(problems,var=[:differences],T=2)

# plot(problems[[2;3]],var=[:difference,:difference_velocity],T=1,legend=:bottomleft)
# plot(problems[[2;6]],var=[:difference,:difference_velocity],T=1,legend=:bottomleft)
# plot(problems[[2;9]],var=[:difference,:difference_velocity],T=1,legend=:bottomleft)
# plot!(problems[[2;4]],var=[:difference,:difference_velocity],T=1,legend=:bottomleft)
# plot!(problems[[2;7]],var=[:difference,:difference_velocity],T=1,legend=:bottomleft)
# plot!(problems[[2;10]],var=[:difference,:difference_velocity],T=1,legend=:bottomleft)
# plot!(problems[[2;5]],var=[:difference,:difference_velocity],T=1,legend=:bottomleft)
# plot(problems[[2;8]],var=[:difference,:difference_velocity],T=1,legend=:bottomleft)
# plot!(problems[[2;11]],var=[:difference,:difference_velocity],T=1,legend=:bottomleft)


# plot!(problems[[1;2]],var=[:difference,:difference_velocity],T=1,legend=:bottomleft)


# plot(problems[[3;4]],var=[:differences])

# plot(problems[[2;3;6]],var=[:differences],T=2,legend=:bottomleft)

# plot(problems[[2;3]],var=[:differences],T=2)
# plot!(problems[[2;4]],var=[:differences],T=2)

# plot(problems[[2;4]],var=[:differences],T=2)
# plot!(problems[[2;5]],var=[:differences],T=2)



# @gif for time in LinRange(0,2,101)
#     plt=plot([pFG3a,pGN],var=[:differences],T=time,legend=:bottomleft,xlims=(0,10),ylims=(-1e-5,1e-5));
#     plot!(plt,[pFG3b,pGN],var=[:differences],T=time,legend=:bottomleft,xlims=(0,10),ylims=(-1e-5,1e-5));

# end

# # Attention au scaling a revoir !
# # Faire la derivee en temps



# function norm(p;T=Inf::AbstractFloat,a,μ,s=0)
#     T=min(max(T,0),p.times.ts[end])
# 	index = findfirst(p.times.ts.>=T)
# 	η,v,u,p,w,x = (p.model.mapfrofull)(p.data.U[index])
# 	mesh=Mesh(x)
#     ∂ₓ=1im*mesh.k
#     |(x)=maximum(abs.(x))
    
#     if s==0
#         return |(η),|(u),|(p),|(w)
#     else
#         return |(ifft(∂ₓ.^s.*fft(η))),|(ifft(∂ₓ.^s.*fft(u))),|(ifft(∂ₓ.^s.*fft(p))),|(ifft(∂ₓ.^s.*fft(w)))
#     end

# end

# function norms(p;T=nothing,a,μ,s=0,rel=false)
#     if T == nothing T = p.times.ts end
#     N=zeros(length(T),4)
#     j=0
#     if rel == true N0 = sum(norm(p;T=T[1],a=a,μ=μ,s=s)) else N0=1 end

#     for t in T
#         j+=1
#         n=norm(p;T=t,a=a,μ=μ,s=s)
#         N[j,:].=n./N0
#     end
#     return N
# end

# function norm(p1,p2;T=Inf::AbstractFloat,a,μ)
#     T=min(max(T,0),p1.times.ts[end],p2.times.ts[end])
# 	i1 = findfirst(p1.times.ts.>=T)
# 	η1,v1,u1 = (p1.model.mapfrofull)(p1.data.U[i1])
# 	i2 = findfirst(p2.times.ts.>=T)
# 	η2,v2,u2 = (p2.model.mapfrofull)(p2.data.U[i2])
# 	mesh=Mesh(x)
#     ∂ₓ=1im*mesh.k
#     |(x)=maximum(abs.(x))
    
#     return |(η1-η2),|(u1-u2)
    
# end

# function norms(p1,p2;T=nothing,a,μ,rel=false)
#     if T == nothing T = p1.times.ts end
#     N=zeros(length(T),2)
#     j=0
#     if rel == true N0 = sum(norm(p1,p2;T=T[1],a=a,μ=μ)) else N0=1 end

#     for t in T
#         j+=1
#         n=norm(p1,p2;T=t,a=a,μ=μ)
#         N[j,:].=n./N0
#     end
#     return N
# end


# plot(norms(pFG3b;a=100,μ=0.01,s=2),label=["η" "u" "p" "w"])

# plot(norms(pGN,pFG1c;a=100,μ=0.01),label=["Δη" "Δu"])


# errid0=[sum(norm(pGN,pFG1a;a=100,μ=0.01)),
#     sum(norm(pGN,pFG2a;a=100,μ=0.01)),
#     sum(norm(pGN,pFG3a;a=100,μ=0.01))]
# errid1=[sum(norm(pGN,pFG1b;a=100,μ=0.01)),
#     sum(norm(pGN,pFG2b;a=100,μ=0.01)),
#     sum(norm(pGN,pFG3b;a=100,μ=0.01))]
# errid2=[sum(norm(pGN,pFG1c;a=100,μ=0.01)),
#     sum(norm(pGN,pFG2c;a=100,μ=0.01)),
#     sum(norm(pGN,pFG3c;a=100,μ=0.01))]


# scatter([20,100,200],errid0,label="no preparation")
# scatter!([20,100,200],errid1,label="w well-prepared")
# scatter!([20,100,200],errid2,label="w and p well-prepared")
# plot!(xscale=:log10,yscale=:log10,legend=:bottomleft,
#         xlabel="a",ylabel="difference (L∞)",
#         title="difference between GN and LCT")


# scatter!([(200,sum(norm(pGN,pFG3c;a=100,μ=0.01)))])
# scatter!([(20,sum(norm(pGN,pFG1c;a=20,μ=0.01)[1]))])
# norm(pGN,pFG2c;a=100,μ=0.01)

# savefig("err.pdf")C