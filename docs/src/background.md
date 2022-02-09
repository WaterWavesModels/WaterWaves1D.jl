# Background

## Water waves

The propagation of waves at the surface of a layer of water is typically modelled
using the incompressible Euler equations inside the fluid domain,
and suitable boundary conditions at the boundaries
(accounting for the impermeable bottom and the free surface).

In an idealized situation where the only external force acting on the fluid is
due to the (constant) vertical gravity acceleration, that the fluid is *homogeneous*
and the flow is *potential* (that is irrotational), the system can be written in closed form
as two evolution equations.

Restricting to *unidimensional waves* (horizontal dimension `d=1`)
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
  ∂_tη-\tfrac{1}{μν} G^μ[ϵη]ψ=0,\\
  ∂_tψ+η+\frac{ϵ}{2ν}|∇ψ|^2-\tfrac{ϵμ}{2ν}\frac{(\tfrac{1}{μν} G^μ[ϵη]ψ+ϵ∇η⋅∇ψ)^2}{1+μϵ²|∇η|^2}=0,
  \end{array}\right.
```
where, by definition,
```math
G^μ[ϵη]ψ=\big(∂_z\Phi-μϵ∇η⋅∇_xΦ\big)\big\vert_{z=ϵη}
```
with ``Φ`` being the unique solution to
```math
\left\{\begin{array}{ll}
μ Δ_x Φ + ∂_z^2 Φ=0& \text{ in } \{(x,z)\ : \  -1<z<ϵη(x)\} , \\
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
[Craig-Sulem](https://doi.org/10.1006/jcph.1993.1164))
is very elegant, it is not directly suitable for efficient numerical simulations,
due to the costly time-dependent elliptic problem involved in the Dirichlet-to-Neumann operator, ``G^μ``.
In our unidimensional framework, it is possible to make use of conformal mapping so as to rewrite the system
using only pointwise operations or Fourier multipliers (that is pointwise operations in Fourier space).
This allows, by means of the [Fast Fourier Transform](https://en.wikipedia.org/wiki/Fast_Fourier_transform) (FFT) and its inverse (IFFT), to provide a very efficient strategy for the numerical simulation of the water waves system. This strategy has been described for instance in [Dyachenko et al.](https://doi.org/10.1016/0375-9601(96)00417-3) and [Choi-Camassa](https://doi.org/10.1061/(asce)0733-9399(1999)125:7(756)).

The resulting code is [`WaterWaves`](@ref WaterWaves1D.WaterWaves).


## Models

Because the above method is relatively recent (in comparison with early studies on water waves),
imperfect (it suffers from "anti-resolution" for large-amplitude waves: the location of gridpoints spread out near wave crests, which in practice may demand the use of a very large number of modes to resolve the flow accurately), and restricted to unidimensional waves, many simplified models have been introduced in the literature. It is the aim of this package to provide a home for many of them.


**Under construction**

### Shallow water models

#### the Saint-Venant system

### Deep water models


## Pseudospectral methods

All the systems described above are discretized (in the space variable) using the Fourier-based pseudospectral method. This method is particularly suitable for data which are either periodic or decaying at infinity, and fairly regular. In this framework one approaches data by a sum [...]
