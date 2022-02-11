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
  ∂_tψ+η+\frac{ϵ}{2ν}|∂_xψ|^2-\tfrac{ϵμ}{2ν}\frac{(\frac{1}{μν} G^μ[ϵη]ψ+ϵ(∂_xη)(∂_xψ))^2}{1+μϵ²|∂_xη|^2}=0,
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
  ∂_tη+∂_x(F₁v + ϵ F₂ (η F₂v))=0,\\[1ex]
  ∂_tv+∂_xη+\tfrac{ϵ}{2}∂_x((F₂v)^2) =0,
  \end{array}\right.
```
with ``F₁=\frac{\tanh(\sqrt\mu D)}{\sqrt\mu D}``, and ``F₂=(F₁)^α`` (here we use the notation ``F(D)`` for the [action](https://en.wikipedia.org/wiki/Multiplier_(Fourier_analysis)) of pointwise multiplying by the function ``F`` in the Fourier space).
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
  ∂_tv+∂_x\big(η+ϵ uv - \tfrac{ϵ}{2}u^2-\tfrac{μϵ}2 (h F₀∂_xu)^2\big) =0,
  \end{array}\right.
```
 where ``h=1 + ϵ η`` is the depth, ``v=∂_xψ`` the derivative of the trace of the velocity potential at the surface, and ``u`` the layer-averaged horizontal velocity obtained by solving the elliptic problem
```math
hu -\tfrac{μ}{3}F₀∂_x( h^3 F₀∂_xu) = hv.
```
with ``F₀=\sqrt{3(F₁^{-1}(D) - 1)}/D`` where ``F₁=\frac{\tanh(\sqrt\mu D)}{\sqrt\mu D}``
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

#### The rectified spectral systems

#### The deep quadratic system

#### The Matsuno system

#### The modified Matsuno system





## Pseudospectral methods

All the systems described above are discretized (in the space variable) using the Fourier-based pseudospectral method. This method is particularly suitable for data which are either periodic or decaying at infinity, and fairly regular. In this framework one approaches data by a sum [...]
