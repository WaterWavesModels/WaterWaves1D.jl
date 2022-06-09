---
title: 'WaterWaves1D.jl: A Julia framework to numerically study the propagation of unidimensional surface gravity waves'
tags:
  - Julia
  - water waves
authors:
  - name: Vincent Duchêne
    orcid: 0000-0002-8349-1284
    equal-contrib: true
    affiliation: 1
  - name: Pierre Navaro
    orcid: 0000-0002-7372-3221
    equal-contrib: true
    affiliation: 1
affiliations:
  - name: Univ Rennes, CNRS, IRMAR - UMR 6625, France.
    index: 1
date: 8 June 2022
bibliography: paper.bib

---

# Summary

The propagation of water waves is a physical phenomenon that is
difficult to calculate by numerical simulations. There are many
models to solve the equations that model this phenomenon. Each model
has advantages and disadvantages. It is interesting to be able to
compare them easily from a single interface.


# Statement of need

WaterWaves1D.jl is a Julia package providing a framework to study
and compare several models for the propagation of unidimensional
surface gravity waves.

Several models are already implemented, including the so-called
water waves system, its truncated spectral expansion, the Green-Naghdi
system, the Matsuno system, and so on.


# Mathematics

The propagation of waves at the surface of a layer of water is
typically modelled using the incompressible Euler equations inside
the fluid domain, and suitable boundary conditions at the boundaries
(accounting for the impermeable bottom and the free surface).

In an idealized situation where the only external force acting on
the fluid is due to the (constant) vertical gravity acceleration,
that the fluid is homogeneous and the flow is potential (that is
irrotational), the system can be written in closed form as two
evolution equations.

$$
  \begin{cases}
  \partial_t\eta-\tfrac{1}{\mu\nu} G^\mu[\epsilon\eta]\psi=0,\\
  \partial_t\psi+\eta+\frac{\epsilon}{2\nu}(\partial_x\psi)^2-\tfrac{\epsilon\mu}{2\nu}\frac{(\frac{1}{\mu} G^\mu[\epsilon\eta]\psi+\epsilon(\partial_x\eta)(\partial_x\psi))^2}{1+\mu\epsilon^2(\partial_x\eta)^2}=0,
  \end{cases}
$$

where, by definition,

$$
G^\mu[\epsilon\eta]\psi=\big(\partial_z\Phi-\mu\epsilon(\partial_x\eta)(\partial_x\phi)\big)\big\vert_{z=\epsilon\eta}
$$

with $\phi$ being the unique solution to the elliptic boundary value problem

$$
\begin{cases}
\mu \partial_x^2 \phi + \partial_z^2 \phi=0& \text{ in } \{ (x,z) \colon   -1<z<\epsilon\eta(x) \} , \\
 \phi = \psi & \text{ on } \{ (x,z) \colon   z=\epsilon\eta(x) \} ,\\
\partial_z \phi =0 & \text{ on } \{ (x,z) \colon z=-1 \} .
\end{cases}
$$

In the above formula,

- $t,x,z$ represent respectively the (rescaled) time, horizontal and vertical space variables.
- $\eta$ represents the *surface deformation*: the free surface may be parametrized as $\{(x,z) :  z=\epsilon\eta(x)\}$.
- $\psi$ is the trace of the *velocity potential* at the surface. In models we generally prefer to use $v=\partial_x\psi$ as the second unknown, since $\psi$ is not necessarily decaying at (spatial) infinity even in finite-energy situations.
- $\epsilon$ is the *nonlinearity* dimensionless parameter, defined as the ratio of the maximal amplitude of the wave to the depth of the layer.
- $\mu$ is the *shallowness* dimensionless parameter, defined as the square of the ratio of the depth of the layer to the typical horizontal wavelength of the flow.
- $\nu$ is a scaling parameter: in shallow water situations one typically sets $\nu=1$ while in deep water situations it is wise to set $\nu=1/\sqrt{\mu}$. In the latter case, $\epsilon\sqrt{\mu}$ being the *steepness* of the wave plays an important role. Especially, taking formally the limit $\mu \rightarrow \infty$ one obtains the infinite-depth situation where the wave steepness is the only remaining parameter.


# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:

- @Duchene:2020  ->  "Vincent Duchêne and Tatsuo Iguchi (2020)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# References
