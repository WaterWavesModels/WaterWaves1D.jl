# Background

```@contents
Depth = 4
Pages = [ "background.md"  ]
```


## Water waves

The propagation of waves at the surface of a layer of water is typically modelled
using the incompressible Euler equations inside the fluid domain,
and suitable boundary conditions at the boundaries
(accounting for the impermeable bottom and the free surface).

In an idealized situation where the only external force acting on the fluid is
due to the (constant) vertical gravity acceleration, that the fluid is *homogeneous*
and the flow is *potential* (that is irrotational), the system can be written in closed form
as two evolution equations.

Restricting to *unidimensional waves* (horizontal dimension ``d=1``)
in a domain without horizontal boundaries,
assuming that the free surface can be described as the graph of a function
which never touches the bottom (that is, non-cavitation)
neglecting surface tension effects,
assuming constant atmospheric pressure at the free surface and
free slip boundary condition at the flat bottom,
the equations in dimensionless variables read as follows
(following notations in [Lannes2013](@citet),
[Duchene2021](@citet)).

```math
  \left\{\begin{array}{l}
  ∂_tη-\tfrac{1}{μν} G^μ[ϵη]ψ=0,\\[1ex]
  ∂_tψ+η+\frac{ϵ}{2ν}(∂_xψ)^2-\tfrac{ϵμ}{2ν}\frac{(\frac{1}{μ} G^μ[ϵη]ψ+ϵ(∂_xη)(∂_xψ))^2}{1+μϵ²(∂_xη)^2}=0,
  \end{array}\right.
```
where, by definition,
```math
G^μ[ϵη]ψ=\big(∂_z\Phi-μϵ(∂_xη)(∂_xΦ)\big)\big\vert_{z=ϵη}
```
with ``Φ`` being the unique solution to the elliptic boundary value problem
```math
\left\{\begin{array}{ll}
μ ∂_x^2 Φ + ∂_z^2 Φ=0& \text{ in } \{(x,z)\ : \  -1<z<ϵη(x)\} , \\
 Φ= ψ & \text{ on } \{(x,z)\ : \  z=ϵη(x)\} ,\\
∂_z Φ=0 & \text{ on } \{(x,z)\ : \  z=-1\} .
\end{array}\right.
```
In the above formula,
* ``t,x,z`` represent respectively the (rescaled) time, horizontal and vertical space variables.
* ``η`` represents the *surface deformation*: the free surface may be parametrized as `` \{(x,z) :  z=ϵη(x)\}``.
* ``ψ`` is the trace of the *velocity potential* at the surface. In models we generally prefer to use ``v=∂ₓψ`` as the second unknown, since ``ψ`` is not necessarily decaying at (spatial) infinity even in finite-energy situations.
* ``ϵ`` is the *nonlinearity* dimensionless parameter, defined as the ratio of the maximal amplitude of the wave to the depth of the layer.
* ``μ`` is the *shallowness* dimensionless parameter, defined as the square of the ratio of the depth of the layer to the typical horizontal wavelength of the flow.
* ``ν`` is a scaling parameter: in shallow water situations one typically sets ``ν=1`` while in deep water situations it is wise to set ``ν=1/\sqrt{μ}``. In the latter case, ``ϵ\sqrt{μ}`` being the *steepness* of the wave plays an important role. Especially, taking formally the limit ``μ→∞`` one obtains the infinite-depth situation where the wave steepness is the only remaining parameter.

While the above formulation (due to [Zakharov1968](@citet)
and [CraigSulem1993](@citet)) is
very elegant, it is not directly suitable for efficient numerical
simulations, due to the costly time-dependent elliptic problem
involved in the Dirichlet-to-Neumann operator, ``G^μ``.  In our
unidimensional framework, it is possible to make use of conformal
mapping so as to rewrite the system using only pointwise operations
or Fourier multipliers (that is pointwise operations in Fourier
space).  This allows, by means of the [Fast Fourier
Transform](https://en.wikipedia.org/wiki/Fast_Fourier_transform)
(FFT) and its inverse (IFFT), to provide a very efficient strategy
for the numerical simulation of the water waves system. This strategy
has been described for instance in [DyachenkoKuznetsovSpectorEtAl1996](@citet) and [ChoiCamassa1999](@citet).

The resulting code is [`WaterWaves`](@ref WaterWaves1D.WaterWaves).

## Models

Because the above method is relatively recent (in comparison with
early studies on water waves), imperfect (it suffers from
"anti-resolution" for large-amplitude waves: the location of
gridpoints spread out near wave crests, which in practice may demand
the use of a very large number of modes to resolve the flow
accurately), and restricted to unidimensional waves, many simplified
models have been introduced in the literature. It is the aim of
this package to provide a home for some of them.


### Shallow water models

Shallow water models are expected to provide valid approximation
to the water waves system for small values of the shallowness
parameter, ``μ≪1`` (in this case, ``ν=1``).

Many of these models are derived and discussed in
[Duchene2021](@citet).

#### The Saint-Venant system

The [Saint-Venant1871](@citet) system
is one of the oldest model for the propagation of water waves. It reads
```math
  \left\{\begin{array}{l}
  ∂_tη+∂_x((1+ϵη)v)=0,\\[1ex]
  ∂_tv+∂_xη+ϵv∂_xv=0.
  \end{array}\right.
```
Notice that the parameter ``μ`` has disappeared:
the Saint-Venant system is the order-zero shallow water model for water waves.

The associated code is [`SaintVenant`](@ref WaterWaves1D.SaintVenant) (see also [`SaintVenant_fast`](@ref WaterWaves1D.SaintVenant_fast)).

#### The Boussinesq systems

The ``abcd``-Boussinesq systems
refers to the full class of equations described in
[BonaChenSaut2002](@citet):

```math
  \left\{\begin{array}{l}
  ∂_tη+∂_x((1+ϵη)u)+a∂_x^3u-b∂_x^2∂_tη=0,\\[1ex]
  ∂_tu+∂_xη+ϵu∂_xu+c∂_x^3η-d∂_x^2∂_tu=0,
  \end{array}\right.
```
where ``a,b,c,d`` can be freely chosen as long as they satisfy ``a+b+c+d=1/3``.

The associated code [`Boussinesq`](@ref WaterWaves1D.Boussinesq)
is restricted to the so-called Hamiltonian case, ``d=b``, and ``c=0``.

Indeed, in that case, one may genuinely interpret ``v=u-d∂_x^2u`` as an approximation to ``∂_xψ``, the derivative of the trace of the velocity potential at the surface, and write the system as evolution equations for the variables ``(η,v)``:
```math
  \left\{\begin{array}{l}
  ∂_tη+∂_x((1+a∂_x^2)(1-b∂_x^2)^{-2}v + ϵ(1-b∂_x^2)^{-1} (η (1-b∂_x^2)^{-1}v))=0,\\[1ex]
  ∂_tv+∂_xη+\tfrac{ϵ}{2}∂_x(((1-b∂_x^2)^{-1}v)^2) =0.
  \end{array}\right.
```


#### The Whitham-Boussinesq systems

The Whitham-Boussinesq systems
can be viewed as modified Boussinesq systems in view of
fully recovering the [dispersive properties](https://en.wikipedia.org/wiki/Dispersion_(water_waves)) of the water waves system. In other words, in the linear framework, that is setting ``ϵ=0``, the model coincides with the linearized water waves system.

Specifically, we consider systems of the form
```math
  \left\{\begin{array}{l}
  ∂_tη+∂_x(F_1^μv + ϵ F_2^μ (η F_2^μv))=0,\\[1ex]
  ∂_tv+∂_xη+\tfrac{ϵ}{2}∂_x((F_2^μv)^2) =0,
  \end{array}\right.
```
with ``F_1^μ=\frac{\tanh(\sqrtμ D)}{\sqrtμ D}``, and ``F_2^μ=(F_1^μ)^α`` (here we use the notation ``F(D)`` for the [action](https://en.wikipedia.org/wiki/Multiplier_(Fourier_analysis)) of pointwise multiplying by the function ``F`` in the Fourier space).
The case ``α = 1`` has been introduced by [DinvayDutykhKalisch2019](@citet), more general situations have been studied by [Emerald2021](@citet).


The associated code is [`WhithamBoussinesq`](@ref WaterWaves1D.WhithamBoussinesq).

#### The Green-Naghdi system

The (Serre-)Green-Naghdi system ([Serre1953](@citet), [SuGardner1969](@citet), [GreenNaghdi1976](@citet))
is sometimes called "fully nonlinear Boussinesq system"
and is expected to provide a better approximation when the parameter ``ϵ`` is large.

One of its many formulations is
```math
  \left\{\begin{array}{l}
  ∂_tη+∂_x\big( h u\big)=0,\\[1ex]
  ∂_tv+∂_x\big(η+ϵ uv - \tfrac{ϵ}{2}u^2-\tfrac{μϵ}2 (h∂_xu)^2\big) =0,
  \end{array}\right.
```
where ``h=1 + ϵ η`` is the depth, ``v=∂_xψ`` the derivative of the trace of the velocity potential at the surface, and ``u`` the layer-averaged horizontal velocity obtained by solving the elliptic problem
```math
 hu -\tfrac{μ}{3}∂_x( h^3 ∂_xu) = hv.
```

The associated code is [`SerreGreenNaghdi`](@ref WaterWaves1D.SerreGreenNaghdi).

#### The square-root depth system

The "√D" system
proposed by [CotterHolmPercival2010](@citet)
can be written as
```math
  \left\{\begin{array}{l}
  ∂_tη+∂_x\big( h u\big)=0,\\[1ex]
  ∂_tv+∂_x\big(η+\tfrac{ϵ}{2}v^2+\tfrac{μϵ}{6h^2} (∂_x(hu))^2\big) =0,
  \end{array}\right.
```
where ``h=1 + ϵ η`` is the depth, ``v=∂_xψ`` the derivative of the trace of the velocity potential at the surface, and ``u`` is obtained by solving the elliptic problem
```math
 u -\tfrac{μ}{3}∂_x( h^{-1} ∂_x(hu)) = v.
```

The associated code is [`SquareRootDepth`](@ref WaterWaves1D.SquareRootDepth).

#### The "non-hydrostatic" system

The "non-hydrostatic" system
proposed by [BristeauMangeneySainte-MarieEtAl2015](@citet)
can be written as
```math
  \left\{\begin{array}{l}
  ∂_tη+∂_x\big( h u\big)=0,\\[1ex]
  ∂_tv+∂_x\big(η+\tfrac{ϵ}{2}v^2+\tfrac{μϵ}{2} ( w∂_x^2(hu) + ∂_x^2( hu w ) + ϵ (∂_xη)(∂_xw)u )\big) =0,
  \end{array}\right.
```
where ``h=1 + ϵ η`` is the depth, ``v=∂_xψ`` the derivative of the trace of the velocity potential at the surface, ``w=-(h ∂ₓu)/2`` and ``u`` is obtained by solving the elliptic problem
```math
 hu -\tfrac{μ}{4}∂_x( h^3 ∂_xu) = hv.
```

The associated code is [`NonHydrostatic`](@ref WaterWaves1D.NonHydrostatic).

#### The Whitham-Green-Naghdi system

The Whitham-Green-Naghdi system
proposed by [DucheneIsrawiTalhouk2015](@citet).
can be viewed as a modified Green-Naghdi system in view of
fully recovering the [dispersive properties](https://en.wikipedia.org/wiki/Dispersion_(water_waves)) of the water waves system. In other words, in the linear framework, that is setting ``ϵ=0``, the model coincides with the linearized water waves system.

```math
  \left\{\begin{array}{l}
  ∂_tη+∂_x\big( h u\big)=0,\\[1ex]
  ∂_tv+∂_x\big(η+ϵ uv - \tfrac{ϵ}{2}u^2-\tfrac{μϵ}2 (h F_0^μ∂_xu)^2\big) =0,
  \end{array}\right.
```
 where ``h=1 + ϵ η`` is the depth, ``v=∂_xψ`` the derivative of the trace of the velocity potential at the surface, and ``u`` the layer-averaged horizontal velocity obtained by solving the elliptic problem
```math
hu -\tfrac{μ}{3}F_0^μ∂_x( h^3 F_0^μ∂_xu) = hv,
```
where ``F_0^μ=\sqrt{3((F_1^μ)^{-1}(D) - 1)}/D`` and ``F_1^μ=\frac{\tanh(\sqrtμ D)}{\sqrtμ D}``
(here we use the notation ``F(D)`` for the [action](https://en.wikipedia.org/wiki/Multiplier_(Fourier_analysis)) of pointwise multiplying by the function ``F`` in the Fourier space).

The associated code is [`WhithamBoussinesq`](@ref WaterWaves1D.WhithamBoussinesq).

#### The Isobe-Kakinuma systems

The Isobe-Kakinuma model is a hierarchy of systems proposed by
[Isobe1994](@citet) and studied in [Kakinuma2001](@citet),
depending on the rank ``N`` and the parameters ``(p_0,p_1,⋯,p_N)``

```math
  \left\{\begin{array}{l}
  ∂_tη+∂_x\left( \sum_{j=0}^N\tfrac{h^{p_j+1}}{p_j+1}∂_xϕ_j \right)=0,\\[1ex]
  ∂_tv+∂_x\left( η
  +ϵ \left( \sum_{i=0}^Np_ih^{p_i-1}ϕ_i \right)∂_x\left( \sum_{j=0}^N\tfrac{h^{p_j+1}}{p_j+1}∂_xϕ_j \right)
  +\tfrac{ϵ}{2}\left( \sum_{j=0}^Nh^{p_j}∂_xϕ_j\right)^2
  +\tfrac{ϵ}{2μ} \left( \sum_{j=0}^Np_jh^{p_j-1}ϕ_j\right)^2 \right) =0,
  \end{array}\right.
```
where ``h=1 + ϵ η`` is the depth, ``v=∂_xψ`` the derivative of the trace of the velocity potential at the surface, and ``(ϕ₀,ϕ₁,⋯,ϕ_N)`` are obtained by solving the elliptic problem
```math
  \left\{\begin{array}{l}
\sum_{j=0}^Nh^{p_j}ϕ_j =ψ,\\[1ex]
-h^{p_i} ∂_x\left(\sum_{j=0}^N\tfrac{h^{p_j+1}}{p_j+1}∂_xϕ_j \right)
+ ∂_x\left(\sum_{j=0}^N\tfrac{h^{p_i+p_j+1}}{p_i+p_j+1}∂_xϕ_j \right)
-\tfrac{1}{μ} \sum_{j=0}^N \tfrac{p_ip_j}{p_i+p_j+1}ϕ_j=0 \quad (\forall i∈\{1,⋯,N\})
  \end{array}\right.
```

The associated code [`IsobeKakinuma`](@ref WaterWaves1D.IsobeKakinuma)
is (for now) limited to the case ``N=1`` and ``(p₀,p₁)=(0,2)``.

#### The Choi systems

The Choi model is a hierarchy of systems proposed by
[Choi2022](@citet),
depending on the rank ``M``. 
We actually consider an asymptotically equivalent variant which reads

```math
  \left\{\begin{array}{l}
∂_tη+∂_x\left(\sum_{m=0}^M  h^{2m+1}\frac{(-μ∂_x^2)^mv }{(2m+1)!}\right) =0,\\[1ex]
\big(	1-\sum_{m=1}^M μ∂_x( h^{2m}∂_x\frac{(-μ∂_x^2)^{m-1}}{(2m)!})\big)∂_tv +∂_xη  \\[1ex] 
		\qquad +∂_x\left(\frac12\sum_{m=0}^M h^{2m}\left(\sum_{j=0}^m \frac{(-μ∂_x^2)^j v}{(2j)!}  \frac{(-μ∂_x^2)^{m-j} v}{(2m-2j)!}-μ\sum_{j=0}^{m-1} \frac{∂_x(-μ∂_x^2)^j v}{(2j+1)!} \frac{∂_x(-μ∂_x^2)^{m-j-1} v}{(2m-2j-1)!} \right)\right),
  \end{array}\right.
```
where ``h=1 + ϵ η`` is the depth and ``v=∂_xϕᵦ`` is the derivative of the trace of the velocity potential at the bottom.

The associated code is [`Choi`](@ref WaterWaves1D.Choi).


### Small steepness models

Small-steepness models rely on the smallness of the steepness dimensionless parameter, ``ϵ\sqrt{μ}≪1``,
and may be valid in shallow water as well as deep water configurations.

In what follows we use the notation ``F(D)`` for the [action](https://en.wikipedia.org/wiki/Multiplier_(Fourier_analysis)) of pointwise multiplying by the function ``F`` in the Fourier space. Specifically,
```math
T^μ=-{\rm i}\tanh(\sqrtμ D)
```
is the "Tilbert transform"
(related to the [Hilbert transform](https://en.wikipedia.org/wiki/Hilbert_transform#Relationship_with_the_Fourier_transform), the latter arising in the infinite layer configuration, ``μ=∞``).

#### The Airy equations

The simplest small-steepness model is the linear [Airy](https://en.wikipedia.org/wiki/Airy_wave_theory) water waves obtained by setting ``ϵ=0`` in the water waves system:
```math
  \left\{\begin{array}{l}
  ∂_tη-\tfrac{1}{\sqrtμ ν} T^μv =0,\\[1ex]
  ∂_tv+∂_xη=0,
  \end{array}\right.
```
where we denote ``v=∂_xψ`` the derivative of the trace of the velocity potential at the surface.

The associated code is [`Airy`](@ref WaterWaves1D.Airy).


#### The spectral systems

[DommermuthYue1987](@citet),
[WestBruecknerJandaEtAl1987](@citet), and
[CraigSulem1993](@citet) have proposed
a hierarchy of systems based on a "spectral" expansion, which can be interpreted
through the Taylor expansion of the Dirichlet-to-Neumann, ``G^μ[ϵη]ψ``, with respect
to the surface deformation variable:
```math
G^μ[ϵη]ψ=G^μ[0]ψ + ϵ (D_η G^μ[0])(ϵη)ψ + ϵ^2 (D_η^2 G^μ[0])(ϵη,ϵη)ψ + ⋯
```

The first nonlinear system of the hierarchy, incorporating only quadratic nonlinearities, is
```math
  \left\{\begin{array}{l}
  ∂_tη-\tfrac{1}{\sqrtμ ν} T^μv  + \tfrac{ϵ}{ν} ∂_x\big(η v +  T^μ(η T^μ v)\big) =0,\\[1ex]
  ∂_tv+∂_xη+\frac{ϵ}{2ν}∂_x\big( v^2-(T^μv)^2\big)=0,
  \end{array}\right.
```
where we denote ``v=∂_xψ`` the derivative of the trace of the velocity potential at the surface.

Higher order systems can be constructed using recursive formula.
Explicit expressions up to quintic nonlinearities are given in
[Choi2019](@citet).

The associated code is [`WWn`](@ref WaterWaves1D.WWn).




#### The rectified spectral systems

It turns out the spectral models above suffer from spurious instabilities; see
[AmbroseBonaNicholls2014](@citet).

[DucheneMelinand2024](@citet) proposed a "rectified" quadratic model:
```math
  \left\{\begin{array}{l}
  ∂_tη-\tfrac{1}{\sqrtμ ν} T^μv  + \tfrac{ϵ}{ν} ∂_x\big((J^δη) v +  T^μ((J^δη) T^μ v)\big) =0,\\[1ex]
  ∂_tv+∂_xη+\frac{ϵ}{2ν}∂_xJ^δ\big( v^2-(T^μv)^2\big)=0,
  \end{array}\right.
```
with ``J^δ=J_0(δD)`` where ``J_0(k)`` approaches ``1`` for small wavenumbers, ``k``, and approaches ``0`` for large wavenumbers; and the parameter ``δ`` can be freely chosen, but is typically of the size of ``\tfrac{ϵ}{ν}``. In the associated code, [`WWn`](@ref WaterWaves1D.WWn), one has by default
```math
J_0(k)=\min(1,1/|k|).
```


#### The Matsuno system

The model introduced by [Matsuno1992](@citet) is
```math
  \left\{\begin{array}{l}
  ∂_tη-\tfrac{1}{\sqrtμ ν} T^μu  + \tfrac{ϵ}{ν} ∂_x(η u) +  \tfrac{ϵ}{ν} T^μ(η ∂_x T^μ u) =0,\\[1ex]
  ∂_tu+\big(1-ϵ\sqrtμ T^μ∂_xη\big)∂_xη+\frac{ϵ}{2ν}∂_x\big( u^2\big)=0,
  \end{array}\right.
```
where ``u=∂_xψ-ϵ\sqrtμ(T^μ∂_xψ)(∂_xη)`` represents the horizontal velocity at the free surface.

The associated codes are [`Matsuno`](@ref WaterWaves1D.Matsuno), and [`Matsuno_fast`](@ref WaterWaves1D.Matsuno_fast) for a less human-readable but more efficient version.


#### The modified Matsuno system


In view of ensuring the stability of the equations,
[DucheneMelinand2024](@citet)
proposed a modified Matsuno system:
```math
  \left\{\begin{array}{l}
  ∂_tη-\tfrac{1}{\sqrtμ ν} T^μu  + \tfrac{ϵ}{ν} ∂_x(η u) +  \tfrac{ϵ}{ν} T^μ(η ∂_x T^μ u) =0,\\[1ex]
  ∂_tu+\exp\big(-ϵ\sqrtμ T^μ∂_xη\big)∂_xη+\frac{ϵ}{2ν}∂_x\big( u^2\big)=0.
  \end{array}\right.
```

The associated codes is [`modifiedMatsuno`](@ref WaterWaves1D.modifiedMatsuno).

#### The Akers-Nicholls system

The model introduced in [AkersNicholls2010](@citet)
(see also [ChengGranero-BelinchonShkollerEtAl2019](@citet))
can be written as
```math
  \left\{\begin{array}{l}
  ∂_tη+∂_x m=0,\\[1ex]
  ∂_tm-\tfrac{1}{\sqrtμ ν} T^μ\big(η+\frac{ϵ}{ν}(L^μ m)^2\big)+\frac{ϵ}{ν}\big(η∂_xη+T^μ(η ∂_x T^μ η)\big)=0,
  \end{array}\right.
```
with notations as above, ``L^μ=\frac{ν\sqrtμ D}{\tanh(\sqrtμ D)}`` and where ``m=-\frac1{\sqrtμ ν} T^μψ  + \frac{ϵ}{ν} \big(η ∂_xψ +  T^μ(η T^μ ∂_xψ)\big)`` represents the vertically integrated horizontal momentum.

The associated codes is [`AkersNicholls`](@ref WaterWaves1D.AkersNicholls), and [`AkersNicholls_fast`](@ref WaterWaves1D.AkersNicholls_fast) for a less human-readable but more efficient version.

### Other models

All the models presented so far consists in two evolution equations for the surface elevation `η` and one velocity variable. Our package is not limited to such structure.

#### Unidirectional models (KdV, BBM, Whitham)

The celebrated [KortewegDeVries1895](@citet) (KdV) model was introduced to model water waves.

```math
∂_tη+∂_x η+\tfrac{3ϵ}{2} η ∂_xη + \tfrac{μ}{6}∂_x^3η=0.
```

The [BenjaminBonaMahony1972](@citet) (BBM) model is a well-known asymptotically equivalent variant.
```math
(1-\tfrac{μ}{6}∂_x^2)∂_tη+∂_x η+\tfrac{3ϵ}{2} η ∂_xη =0.
```

A full dispersion model was proposed by [Whitham1967](@citet)
```math
∂_tη+M^μ∂_x η+\tfrac{3ϵ}{2} η ∂_xη =0.
```
where ``M^μ=\big(\tfrac{\tanh(\sqrtμ D)}{\sqrtμ D}\big)^{1/2}``.

While these models have been designed to approximate *unidirectional* waves (right-going or left-going after sign-modifications), they can also be used as building blocks to approximate general water waves; see [Lannes2013](@citet) and [Emerald2021b](@citet).

The associated codes are [`KdV`](@ref WaterWaves1D.KdV), [`BBM`](@ref WaterWaves1D.BBM) and [`Whitham`](@ref WaterWaves1D.Whitham) respectively.

#### Relaxed Green-Naghdi models

[FavrieGavrilyuk2017](@citet), [EscalanteDumbserCastro2019](@citet) and [Richard2021](@citet) among others have proposed a strategy to efficiently approximate solutions of the Green-Naghdi model.

It consists in solving an augmented system with additional unknowns and a relaxation parameter, for instance

```math
  \left\{\begin{array}{l}
  ∂_tη+∂_x (hu)=0,\\[1ex]
  h(∂_tu+ϵu∂_x u+∂_x η)+aμ ∂_x (hp)=0,\\[1ex]
  h(∂_tp +ϵu∂_x p)+a(2w+h∂_xu)=0,\\[1ex]
	h(∂_tw+uϵ∂_x w) = a\tfrac{3}{2}p.
  \end{array}\right.
```
where ``h=1 + ϵ η`` is the water depth and ``p`` and ``w`` are expected to approximate the layer-veraged pressure and vertical velocity if ``a\gg1``.

Notice the system has four evolution equations, in particular initial data for the augmented variables ``p`` and ``w`` must be suitably chosen. 

The associated code is [`relaxedGreenNaghdi`](@ref WaterWaves1D.relaxedGreenNaghdi).

#### Two-dimensional models

Despite its name, our package is not strictly restricted to one-dimensional evolution equations, although the full water waves equations with two horizontal dimensions is largely outside its scope.

As a demonstration we implemented the two-dimensional [Saint-Venant1871](@citet) system
which reads
```math
  \left\{\begin{array}{l}
  ∂_tη+\nabla\cdot((1+ϵη)\bm{v})=0,\\[1ex]
  ∂_tv+\nabla η+ϵ(\bm{v}\cdot\nabla)\bm{v}=0.
  \end{array}\right.
```
Alternatively, we may replace ``(\bm{v}\cdot\nabla)\bm{v}`` with ``\tfrac12 \nabla(|\bm{v}|^2)`` so a to preserve the hamiltonian structure of the original model, and we implemented both versions.


The associated code is [`SaintVenant2D`](@ref WaterWaves1D.SaintVenant2D) (see also [`SaintVenant2D_fast`](@ref WaterWaves1D.SaintVenant2D_fast)).


## Mass, momentum, energy

Solutions to the water waves equations preserve along time (among other integrals of motion listed by [BenjaminOlver1982](@citet))
- the excess of mass,
```math
\int η(t,x)\, {\rm d}x\ ;
```
- the momentum (or rather the horizontal impulse),
```math
\int (η ∂_xψ)(t,x)\, {\rm d}x\ ;
```
- the total energy,
```math
\frac12 \int (η^2 + \tfrac{1}{μν} ψ G^μ[ϵη]ψ)(t,x)\, {\rm d}x\ .
```
Such is also the case for solutions to all models cited above using as velocity variable ``v=∂_xψ`` (with analogous formula).

The numerical discretization of these preserved quantities, coded in [`mass`](@ref WaterWaves1D.mass), [`momentum`](@ref WaterWaves1D.momentum), [`energy`](@ref WaterWaves1D.energy) (and [`mass_diff`](@ref WaterWaves1D.mass_diff), [`momentum_diff`](@ref WaterWaves1D.momentum_diff), [`energy_diff`](@ref WaterWaves1D.energy_diff)  for the difference between initial and final time) provide valuable insights at the precision of a computed numerical solution.




## Pseudospectral methods

Although this is not imperative of the package, all the codes mentioned above use Fourier-based pseudospectral methods for spatial discretization. This method is particularly suitable for data which are either periodic or decaying at infinity (in which case the function at stake is considered as periodic on a sufficiently large period), and fairly regular. In this framework one approaches ``2L``-periodic functions by *finite* Fourier sums of the form (assuming ``N`` even)
```math
u(t,x)≈\sum_{k=-N/2}^{N/2-1} a_k(t) e^{{\rm i} \tfrac{π}{L}k x},
```
and seek a system of ordinary differential equations on the coefficients ``a_k(t)`` (which are then discretized in time by your favorite time solver).

In practice, given the values ``(u(x_j))_{j∈\{0,⋯,N-1\}}`` at the regularly spaced *collocation points* ``x_j= -L+2j\tfrac{L}{N}``, one uses the discrete Fourier transform (computed efficiently with a [Fast Fourier transform](https://en.wikipedia.org/wiki/Fast_Fourier_transform) (FFT)) to deduce the corresponding coefficients ``a_k``. Hence the coefficients are related to the coefficients of the (infinite) Fourier series
```math
u(t,x)=\sum_{k=-∞}^{∞} c_k(t) e^{{\rm i} \tfrac{π}{L}k x}
```
through the relation
```math
a_k(t) = \sum_{ℓ\in\mathbb{Z}} c_{k+ℓN}(t).
```
Hence each of the coefficients ``a_k`` encompasses a full series of Fourier coefficients,
which is unavoidable since ``e^{{\rm i} \tfrac{π}{L}k x}`` and ``e^{{\rm i} \tfrac{π}{L}(k+ℓN) x}`` are indistinguishable on the gridpoints ``x∈\{x_0,⋯,x_{N-1}\}``: this is called [*aliasing*](https://en.wikipedia.org/wiki/Aliasing).
Yet for smooth functions, Fourier coefficients rapidly decrease, and the error between the
finite spectral decomposition and the infinite Fourier decomposition can often be made immaterial (that is to the order of machine precision rounding errors) when choosing a sufficient large number of modes, ``N``.

When using the finite spectral decomposition, the action of [Fourier multipliers](https://en.wikipedia.org/wiki/Multiplier_(Fourier_analysis)), and in particular spatial differentiation, is performed *exactly*
(that is up to machine precision rounding errors) through corresponding multiplication on the discrete coefficients ``a_k``, and only nonlinear contributions require some attention. The simplest way to approximately compute products (or any pointwise operations in the space variable) is by performing
pointwise operations on values at collocation points which are obtained from the discrete coefficients using discrete inverse Fourier transform. For the sake of discussion, consider

```math
 u^2(x_j,t) ≈ \big(\sum_{k=-N/2}^{N/2-1} a_k(t) e^{{\rm i} \tfrac{π}{L}k x_j}\big)^2
 = \sum_{m=-N/2}^{N/2-1}\sum_{n=-N/2}^{N/2-1} a_m a_n e^{{\rm i} \tfrac{π}{L}(m+n) x_j}.
```

One infers

```math
u^2(x_j,t) ≈ \sum_{k=-N/2}^{N/2-1} b_k(t) e^{{\rm i} \tfrac{π}{L}k x_j}, \qquad b_k = \sum_{m+n\in \{k-N,k,k+N\}}a_m a_n.
```

In the above formula for ``b_k``, some of the summands are spurious
effects from aliasing, which sometimes contribute to numerical
instabilities. In order to suppress such terms the so-called
*dealiasing* consists in adding a sufficient number of modes with
coefficients set to zero (in practice one often uses ideal low-pass
filters, that is set to zero extreme modes, so as to always work
with vectors with a fixed given length). The so-called
[Orszag1971](@citet)'s `3/2` rule states that, in the presence of
quadratic nonlinearities, padding `3/2` modes (or zero-ing `1/3`
modes) is sufficient to discard all spurious aliasing contributions,
and provides in particular a
[Galerkin](https://en.wikipedia.org/wiki/Galerkin_method) approximation
since the error is orthogonal to all expansion functions.

