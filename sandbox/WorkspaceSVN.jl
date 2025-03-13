# # Workspace
#
#md # [`notebook`](@__NBVIEWER_ROOT_URL__notebooks/two_problems.ipynb)
#
using WaterWaves1D, Plots, LinearAlgebra;

gr();

#---- parameters


para = (
    ϵ = 1/4,
    κ = 0.01,
    ν = 0*0.01,
    N = 2^7, # number of collocation points
    L = π,# size of the mesh (-L,L)
    T = 1,# final time of computation
    dt = 0.001,  # timestep
);

hcont(ϱ)=0.5 .+ 8*(ϱ.-1.5).^2
ucont(ϱ)=ϱ-1.5

function phys(M;centered=true) 
    if centered == false
        ρ =  Vector(1 .+ (0:M-1)/M)
    else
        ρ =  Vector(1 .+ (0:M)/M)
        ρ = (ρ[1:end-1]+ρ[2:end])/2
    end
    return (
    ρ =  ρ,
    δ = hcont.(ρ),
    u₀ = ucont.(ρ)
        )
end

param(M;centered=true) = merge(para,phys(M;centered=centered))



#---- initial data
function init(M)
    ζM(x) = sin.(x)*(phys(M).δ)';
    uM(x) = zero(x)*(phys(M).δ)';
    return Init(ζM,uM)
end

pb(M;centered=true) = Problem( SVN(param(M;centered=centered); dealias = 1) , init(M), param(M;centered=centered), solver=RK4(param(M),2*M), label = "SV$M")

p4=pb(105)
p3=pb(35)
p2=pb(21)
p1=pb(15)

solve!(p1)
solve!(p2)
solve!(p3)
solve!(p4)


#---- visualization
x=Mesh(para).x
sol1=solution(p1)
sol2=solution(p2)
sol3=solution(p3)
sol4=solution(p4)

err = [
norm(sol3[1][1]-sol4[1][2],Inf);
norm(sol2[1][1]-sol4[1][3],Inf);
norm(sol1[1][1]-sol4[1][4],Inf) ]

N= [35 ; 21; 15]

scatter!(N,err,xaxis=:log10,yaxis=:log10, label = "ρ centré")
plot!(N, 1 ./ N.^2, label = "N^{-2}")


#---- Decentre

p4b=pb(120;centered=false)
p3b=pb(60;centered=false)
p2b=pb(30;centered=false)
p1b=pb(15;centered=false)

solve!(p1b)
solve!(p2b)
solve!(p3b)
solve!(p4b)


sol1b=solution(p1b)
sol2b=solution(p2b)
sol3b=solution(p3b)
sol4b=solution(p4b)

errb = [
norm(sol3b[1][2]-sol4b[1][3],Inf);
norm(sol2b[1][2]-sol4b[1][5],Inf);
norm(sol1b[1][2]-sol4b[1][9],Inf) ]

Nb= [60 ; 30; 15]

scatter(Nb,errb,xaxis=:log10,yaxis=:log10, label = "ρ décentré")
plot!(Nb, 1 ./ Nb, label = "N^{-1}")


h12 = (phys(15).δ)' .+ para.ϵ * sol1[1]
h24 = (phys(21).δ)' .+ para.ϵ * sol2[1]

function H(problem,para;T=Inf) 
    sol = solution(problem;T=T)
    M = size(sol[1])[2]
    h = (phys(M).δ)' .+ para.ϵ * sol[1]
    return cumsum(h[:,end:-1:1],dims=2)[:,end:-1:1]/size(h)[2]
end

plot(x,H(p1,para) ,ylims=(0,1.5),color=:1,label="")
plot!(x,H(p2,para),ylims=(0,1.5),color=:2,label="")



@gif for time in LinRange(0,para.T,101)
    plot(x,H(p1,para;T=time),ylims=(0,1.5),label="")
end
