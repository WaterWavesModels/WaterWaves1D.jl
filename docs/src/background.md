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
(following notations in [Lannes](https://bookstore.ams.org/surv-188),
[Duchêne](https://www.ams.org/open-math-notes/omn-view-listing?listingId=111309)).

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

While the above formulation (due to
[Zakharov](https://doi.org/10.1007/BF00913182) and
[Craig and Sulem](https://doi.org/10.1006/jcph.1993.1164))
is very elegant, it is not directly suitable for efficient numerical simulations,
due to the costly time-dependent elliptic problem involved in the Dirichlet-to-Neumann operator, ``G^μ``.
In our unidimensional framework, it is possible to make use of conformal mapping so as to rewrite the system
using only pointwise operations or Fourier multipliers (that is pointwise operations in Fourier space).
This allows, by means of the [Fast Fourier Transform](https://en.wikipedia.org/wiki/Fast_Fourier_transform) (FFT) and its inverse (IFFT), to provide a very efficient strategy for the numerical simulation of the water waves system. This strategy has been described for instance in [Dyachenko et al.](https://doi.org/10.1016/0375-9601(96)00417-3) and [Choi and Camassa](https://doi.org/10.1061/(asce)0733-9399(1999)125:7(756)).

The resulting code is [`WaterWaves`](@ref WaterWaves1D.WaterWaves).

## Models

Because the above method is relatively recent (in comparison with early studies on water waves),
imperfect (it suffers from "anti-resolution" for large-amplitude waves: the location of gridpoints spread out near wave crests, which in practice may demand the use of a very large number of modes to resolve the flow accurately), and restricted to unidimensional waves, many simplified models have been introduced in the literature. It is the aim of this package to provide a home for some of them.


### Shallow water models

Shallow water models are expected to provide valid approximation to the water waves system
for small values of the shallowness parameter, ``μ≪1`` (in this case, ``ν=1``).

Many of these models are derived and discussed in
[Duchêne](https://www.ams.org/open-math-notes/omn-view-listing?listingId=111309).

#### The Saint-Venant system

The [Saint-Venant](https://en.wikipedia.org/wiki/Shallow_water_equations#One-dimensional_Saint-Venant_equations) system
is one of the oldest model for the propagation of water waves. It reads
```math
  \left\{\begin{array}{l}
  ∂_tη+∂_x((1+ϵη)v)=0,\\[1ex]
  ∂_tv+∂_xη+ϵv∂_xv=0.
  \end{array}\right.
```
Notice that the parameter ``μ`` has disappeared:
the Saint-Venant system is the order-zero shallow water model for water waves.

The associated code is [`SaintVenant`](@ref WaterWaves1D.SaintVenant).

#### The Boussinesq systems

The ``abcd``-[Boussinesq](https://en.wikipedia.org/wiki/Boussinesq_approximation_(water_waves)) systems
refers to the full class of equations described by
[Bona, Chen and Saut](https://doi.org/10.1007/s00332-002-0466-4)

```math
  \left\{\begin{array}{l}
  ∂_tη+∂_x((1+ϵη)u)+a∂_x^3u-b∂_x^2∂_tη=0,\\[1ex]
  ∂_tu+∂_xη+ϵu∂_xu+c∂_x^3η-d∂_x^2∂_tu=0,
  \end{array}\right.
```
where ``a,b,c,d`` can be freely chosen as long as they satisfy ``a+b+c+d=1/3``.

The associated code [`Boussinesq`](@ref WaterWaves1D.Boussinesq)
is restricted to the so-called Hamiltonian case, ``b=d``, and ``c=0``.

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
with ``F_1^μ=\frac{\tanh(\sqrt\mu D)}{\sqrt\mu D}``, and ``F_2^μ=(F_1^μ)^α`` (here we use the notation ``F(D)`` for the [action](https://en.wikipedia.org/wiki/Multiplier_(Fourier_analysis)) of pointwise multiplying by the function ``F`` in the Fourier space).
The case ``α = 1`` has been introduced by [Dinvay, Dutykh and Kalisch](https://doi.org/10.1016/j.apnum.2018.09.016), more general situations have been studied by [Emerald](https://doi.org/10.1137/20M1332049).


The associated code is [`WhithamBoussinesq`](@ref WaterWaves1D.WhithamBoussinesq).

#### The Green-Naghdi system

The (Serre-)Green-Naghdi system ([Serre](https://10.1051/lhb/1953058), [Su and Gardner](https://10.1063/1.1664873), [Green and Naghdi](https://10.1017/s0022112076002425))
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
proposed by [Cotter, Holm and Percival](https://doi.org/10.1098/rspa.2010.0124)
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
proposed by [Bristeau, Mangeney, Sainte-Marie and Seguin](https://doi.org/10.3934/dcdsb.2015.20.961)
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
proposed by [Duchêne, Israwi and Talhouk](https://doi.org/10.1137/130947064).
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
hu -\tfrac{μ}{3}F_0^μ∂_x( h^3 F_0^μ∂_xu) = hv.
```
with ``F_0^μ=\sqrt{3((F_1^μ)^{-1}(D) - 1)}/D`` where ``F_1^μ=\frac{\tanh(\sqrt\mu D)}{\sqrt\mu D}``
(here we use the notation ``F(D)`` for the [action](https://en.wikipedia.org/wiki/Multiplier_(Fourier_analysis)) of pointwise multiplying by the function ``F`` in the Fourier space).

The associated code is [`WhithamBoussinesq`](@ref WaterWaves1D.WhithamBoussinesq).

#### The Isobe-Kakinuma systems

The Isobe-Kakinuma model is a hierarchy of models proposed by
[Isobe](https://doi.org/10.1061/9780784400890.023),
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



### Small steepness models

Small-steepness models rely on the smallness of the steepness dimensionless parameter, ``ϵ\sqrt{μ}≪1``,
and may be valid in shallow water as well as deep water configurations.

**Descriptions under completion**


#### The spectral systems

[Dommermuth and Yue](https://doi.org/10.1017/s002211208700288x),
[West et al.](https://doi.org/10.1029/jc092ic11p11803), and
[Craig and Sulem](https://doi.org/10.1006/jcph.1993.1164) have proposed
a hierarchy of systems based on a "spectral" expansion, which can be interpreted
through the Taylor expansion of the Dirichlet-to-Neumann, ``G^μ[ϵη]ψ``, with respect
to the surface deformation variable:
```math
G^μ[ϵη]ψ=G^μ[0]ψ + ϵ (D_η G^μ[0])(ϵη)ψ + ϵ^2 (D_η^2 G^μ[0])(ϵη,ϵη)ψ + ⋯
```

The first nonlinear system of the hierarchy, incorporating only quadratic nonlinearities, is
```math
  \left\{\begin{array}{l}
  ∂_tη-\tfrac{1}{\sqrt μ ν} T^μv  + \tfrac{ϵ}{ν} ∂_x\big(η v +  T^μ(η T^μ v)\big) =0,\\[1ex]
  ∂_tv+∂_xη+\frac{ϵ}{2ν}∂_x\big( v^2-(T^μv)^2\big)=0,
  \end{array}\right.
```
with ``v=∂_xψ`` the derivative of the trace of the velocity potential at the surface,
and ``T^μ=-i\tanh(\sqrt\mu D)`` the "Tilbert transform"
(related to the [Hilbert transform](https://en.wikipedia.org/wiki/Hilbert_transform#Relationship_with_the_Fourier_transform);
here we use the notation ``F(D)`` for the [action](https://en.wikipedia.org/wiki/Multiplier_(Fourier_analysis)) of pointwise multiplying by the function ``F`` in the Fourier space).

Higher order systems can be constructed using recursive formula.
Explicit expressions up to quintic nonlinearities are given in
[Choi](http://hdl.handle.net/2433/251940).

The associated code is [`WWn`](@ref WaterWaves1D.WWn).




#### The rectified spectral systems

It turns out the spectral models above suffer from spurious instabilities; see
[Ambrose, Bona and Nicholls](https://doi.org/10.1098/rspa.2013.0849).

[Duchêne and Melinand](to be submitted) proposed a "rectified" quadratic model:
```math
  \left\{\begin{array}{l}
  ∂_tη-\tfrac{1}{\sqrt μ ν} T^μv  + \tfrac{ϵ}{ν} ∂_x\big((J^δη) v +  T^μ((J^δη) T^μ v)\big) =0,\\[1ex]
  ∂_tv+∂_xη+\frac{ϵ}{2ν}∂_xJ^δ\big( v^2-(T^μv)^2\big)=0,
  \end{array}\right.
```
with notations as above and ``J^δ=J_0(δD)`` where ``J_0(k)`` approaches ``1`` for small wavenumbers, ``k``, and approaches ``0`` for large wavenumbers; and the parameter ``δ`` can be freely chosen, but is typically of the size of ``\tfrac{ϵ}{ν}``. In the associated code, [`WWn`](@ref WaterWaves1D.WWn), one has by default
```math
J_0(k)=\min(1,1/|k|).
```

#### The deep quadratic system

***Under redaction***

#### The Matsuno system

The model introduced by [Matsuno](https://doi.org/10.1103/PhysRevLett.69.609) is
```math
  \left\{\begin{array}{l}
  ∂_tη-\tfrac{1}{\sqrt μ ν} T^μu  + \tfrac{ϵ}{ν} ∂_x(η u) +  \tfrac{ϵ}{ν} T^μ(η ∂_x T^μ u) =0,\\[1ex]
  ∂_tu+\big(1-ϵ\sqrt μ T^μ∂_xη\big)∂_xη+\frac{ϵ}{2ν}∂_x\big( u^2\big)=0,
  \end{array}\right.
```
where ``u=∂_xψ-ϵ\sqrt μ(T^μ∂_xψ)(∂_xη)`` represents the horizontal velocity at the free surface.

The associated codes are [`Matsuno`](@ref WaterWaves1D.Matsuno), and [`Matsuno_fast`](@ref WaterWaves1D.Matsuno_fast) for a less human-readable but more efficient version.


#### The modified Matsuno system


In view of ensuring the stability of the model,
[Duchêne and Melinand](to be submitted)
proposed a modified Matsuno system:
```math
  \left\{\begin{array}{l}
  ∂_tη-\tfrac{1}{\sqrt μ ν} T^μu  + \tfrac{ϵ}{ν} ∂_x(η u) +  \tfrac{ϵ}{ν} T^μ(η ∂_x T^μ u) =0,\\[1ex]
  ∂_tu+\exp\big(-ϵ\sqrt μ T^μ∂_xη\big)∂_xη+\frac{ϵ}{2ν}∂_x\big( u^2\big)=0.
  \end{array}\right.
```

The associated codes is [`modifiedMatsuno`](@ref WaterWaves1D.modifiedMatsuno).



## Pseudospectral methods

Although this is not imperative of the package, all the codes mentioned above use Fourier-based pseudospectral methods for spatial discretization. This method is particularly suitable for data which are either periodic or decaying at infinity (in which case the function at stake is considered as periodic on a sufficiently large period), and fairly regular. In this framework one approaches ``2L``-periodic functions by *finite* Fourier sums of the form
```math
u(t,x)≈\sum_{k=1}^{N} a_k(t) e^{{\rm i} \tfrac{π}{L}k x},
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
 u^2(x_j,t) ≈ \big(\sum_{k=1}^{N} a_k(t) e^{{\rm i} \tfrac{π}{L}k x_j}\big)^2 = \sum_{m=1}^N\sum_{n=1}^N a_m a_n e^{{\rm i} \tfrac{π}{L}(m+n) x_j}.
```
One infers
```math
u^2(x_j,t) ≈ \sum_{k=1}^{N} b_k(t) e^{{\rm i} \tfrac{π}{L}k x_j}, \qquad b_k = \sum_{m+n-k\in N\mathbb{Z} }a_m a_n.
```
In the above formula for ``b_k``, some of the summands are spurious effects from aliasing, which sometimes contribute to numerical instabilities. In order to suppress such terms the so-called *dealiasing* consists in adding a sufficient number of modes with coefficients set to zero (in practice one often uses ideal low-pass filters, that is set to zero extreme modes, so as to always work with vectors with a fixed given length). The so-called [Orszag](https://doi.org/10.1175/1520-0469(1971)028<1074:oteoai>2.0.co;2)'s `3/2` rule states that, in the presence of quadratic nonlinearities, padding `3/2` modes (or zero-ing `1/3` modes) is sufficient to discard all spurious aliasing contributions, and provides in particular a [Galerkin](https://en.wikipedia.org/wiki/Galerkin_method) approximation since the error is orthogonal to all expansion functions.
