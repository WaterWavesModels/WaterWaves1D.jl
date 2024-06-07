# # Workspace
#
using WaterWaves1D, FFTW;
using Plots; gr();
using DelimitedFiles


#### COMPUTATIONS JULIA

# Construction de la donnee initiale
# Donnees du code d'Arnaud
phys = (Xmin	=	0,
        Xmax	=	200,
        X0		=   100,
        Final_Time	=	10,
        h0      =   1,
        g = 9.81,
        a = 100)

# Parametres de l'onde solitaire de GN
a	    = 0.2* phys.h0
c0	    = sqrt(phys.g * ( phys.h0 + a ))
kappa   = sqrt(3*a/(4*(phys.h0^2)*(phys.h0+a)))

ksi(x)		 = a/(cosh(kappa*(x-phys.X0))^2)
u0(x)		 = c0*(1-phys.h0/(ksi(x)+phys.h0))
h0(x)        = phys.h0+ksi(x)

L=1/kappa#(phys.Xmax-phys.Xmin)
H = phys.h0
c₀=sqrt(phys.g * H)

# Parametres du schema numerique
param = (
        ϵ = 1,
        μ = (H/L)^2,
        a = phys.a/c₀*L/H,
        N = 2^9, # number of collocation points
        xmin = phys.Xmin/L, # size of the mesh (-L,L)
        xmax = phys.Xmax/L,
        T = phys.Final_Time*c₀/L, # final time of computation
        dt = phys.Final_Time*c₀/L*1e-4-eps(),  # timestep
        );

# Construction de la donnee initiale
# # Same thing as 
# (η₁,u₁,v₁,mesh) = SolitaryWaveSerreGreenNaghdi(merge(param,(c=c0/c₀,)); x₀=phys.X0/L)
# init = Init(mesh,η₁,v₁)

mesh=Mesh(param)
k = mesh.k
Dx       =  1im * k
DxF(v) = real.(ifft(Dx .* fft(v)))
η₀ = ksi.(L*mesh.x)/H
h₀ = 1 .+ η₀
u₀ = u0.(L*mesh.x)/c₀
v₀ = u₀ - param.μ/3 ./h₀ .* (DxF(h₀.^3 .*DxF(u₀)))
init = Init(mesh,η₀,v₀) 


#pWW = Problem( WaterWaves(param; dealias = 1) , init, param) 
#pGN = Problem( WhithamGreenNaghdi(param; SGN=true, dealias = 1) , init, param) 
pLCT = Problem( relaxedGreenNaghdi(param; FG=false, id = 0, dealias = 1) , init, param) 

solve!(pLCT);
plot(pLCT) 

# Si besoin, comparaison avec une autre simulation numerique
pFG = Problem( relaxedGreenNaghdi(param; FG=true, id = 0, dealias = 1) , init, param) 

η2,v2,x2=solution(pFG,T=param.T)
η3,v3,x3=solution(pLCT,T=param.T)

eta2 = interpolate(Mesh(x2),η2,x2)
eta3 = interpolate(Mesh(x3),η3,x2)

plot(x2/L,[eta2-eta3])
plot([pLCT2 pLCT3],var=:difference)




#### COMPUTATIONS ARNAUD

# Recupere les donnees du code d'Arnaud

aLCT1 =readdlm("/home/vduchene/LcT_DG_1D/Results_1d/test2_LCT1_40.txt")
aLCT2 =readdlm("/home/vduchene/LcT_DG_1D/Results_1d/test2_LCT2_40.txt")
aLCT4 =readdlm("/home/vduchene/LcT_DG_1D/Results_1d/test2_LCT4_40.txt")
aLCT8 =readdlm("/home/vduchene/LcT_DG_1D/Results_1d/test2_LCT8_40.txt")
#aLCT16 =readdlm("/home/vduchene/LcT_DG_1D/Results_1d/test2_LCT16_4.txt")

scale = 1

Nn1 = Int(length(aLCT1[:,1])/scale)
indices=range(1,step=scale,stop=Nn1)
xx1  =aLCT1[indices,1]
etaLCT1=aLCT1[indices,2].-1.

Nn2 = Int(length(aLCT2[:,1])/scale)
indices=range(1,step=scale,stop=Nn2)
xx2  =aLCT2[indices,1]
etaLCT2=aLCT2[indices,2].-1.

Nn4 = Int(length(aLCT4[:,1])/scale)
indices=range(1,step=scale,stop=Nn4)
xx4  =aLCT4[indices,1]
etaLCT4=aLCT4[indices,2].-1.

Nn8 = Int(length(aLCT8[:,1])/scale)
indices=range(1,step=scale,stop=Nn8)
xx8  =aLCT8[indices,1]
etaLCT8=aLCT8[indices,2].-1.


#Nn16 = Int(length(aLCT16[:,1])/scale)
#indices=range(1,step=scale,stop=Nn16)
#xx16  =aLCT16[indices,1]
#etaLCT16=aLCT16[indices,2].-1.




######## COMPARAISON

η,v,x,=solution(pLCT,T=param.T)
mesh=Mesh(x)
#η2 =   real.(ifft(pLCT.data.U[10001][:,1]))
eta1 = interpolate(mesh,η,xx1/L)
eta2 = interpolate(mesh,η,xx2/L)
eta4 = interpolate(mesh,η,xx4/L)
eta8 = interpolate(mesh,η,xx8/L)
#eta16 = interpolate(mesh,η,xx16/L)



plot( xx1/L,[ etaLCT1/H ],label="LCT1_Arnaud")
plot!(xx2/L,[ etaLCT2/H ],label="LCT2_Arnaud")
plot!(xx4/L,[ etaLCT4/H ],label="LCT4_Arnaud")
plot!(xx8/L,[ etaLCT8/H ],label="LCT8_Arnaud")
plot!(xx8/L, eta8, label = "LCT_Vincent")
ylims!(-0.002,0.001)
#xlims!(0.65,0.75)

plot( xx1/L,[etaLCT1/H - eta1])
plot( xx2/L,[etaLCT2/H - eta2])
plot!(xx4/L,[etaLCT4/H - eta4])
plot!(xx8/L,[etaLCT8/H - eta8])
plot( xx8/L,[ksi.(xx8) - eta8])



plot( xx8/L,[ksi.(xx8) - init.η(xx8/L)])
#plot!(xx8/L,[eta8 - eta8b])
#plot!(xx16/L,[etaLCT16/H - eta16])
ylims!(-1e-4,1e-4)
xlims!(0.4,0.6)

