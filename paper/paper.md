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
date: 03 July 2026
bibliography: paper.bib

---

# Summary

`WaterWaves1D.jl` is a software that aims at facilitating the study and comparison of models for the propagation of water waves in an idealised framework. 
Several models are already implemented, including the so-called water waves system as well as many standard or more advanced simplified models.

`WaterWaves1D.jl` leverages the Julia programming language [@BezansonEdelmanKarpinskiEtAl17] capabilities to offer highly efficient codes with a friendly user interface. The architecture of the package has been designed to ease the inclusion of new models or numerical schemes. 
The package provides useful features for the monitoring and comparison of different models or numerical strategies on a common ground, as well as pleasant illustration options.


# Statement of need

The propagation of waves at the surface of water under the action of gravitational forces (herein referred to as **water waves**) is arguably the physical phenomenon for which the most models have been proposed and are currently in use [@Lannes; @MM4WW]. This is partly due to the variety of regimes describing physical situations —such as shallow water, small-slope, long waves, wave packets— that can be considered and drive the modelling procedure. Even within the same physical regime, different approaches —Boussinesq–Rayleigh expansions, variational methods, spectral methods to name a few— yield different models, and *a posteriori* adjustments such as the BBM trick or the Nwogu approach further expand the range of models proposed in the literature. Yet all these models are built from the same "master equations", and in principle can easily be compared.
Contrarily to softwares that typically focus on a single model or a handful of variations as well as a dedicated numerical strategy, 
the Julia package `WaterWaves1D.jl` provides a unified groundwork for the study and comparison of models.


# Software design

`WaterWaves1D.jl` provides the architecture to solve initial-value problems for PDEs through the method of lines: the full space-time discretization is decomposed into (i) spatial semidiscretization; and (ii) numerical integration of the resulting ODE system. 
Within this framework a numerical problem is defined by the spatial semidiscretization associated with the system of partial differential equations at stake, the initial data and the time integration scheme. In addition, functions mapping back and forth physical quantities to the objects that are numerically computed allow for the comparison and illustration of the numerical results.

A common feature of all the implemented models is that they involve combinations of pointwise operations in physical space, and Fourier multipliers (that is multiplication in Fourier space). In a periodic framework —which can approximate the full real-line framework if the period is sufficiently large— this calls for the use of Fourier spectral methods [@Trefethen00] for spatial semidiscretization. Indeed, the Fast Fourier Transform (FFT) algorithm (and its inverse) allows to convert very efficiently a sequence of values at collocation points to a representation in the frequency domain (and back) and benefits from the spectral (or exponential) convergence property on smooth functions. Yet this is not a restriction of the architecture of the package, and other methods for spatial discretization can be used.

The architecture of the package allows in principle any solver for time integration provided that the semi-discretized model supplies the necessary operators. All implemented models can be equipped with any explicit Runge–Kutta methods, and the standard explicit Euler and fourth order Runge-Kutta schemes are provided. Additionnally, as a proof-of-concept, the symplectic Euler and Störmer–Verlet [@HairerLubichWanner03] and the exponential Euler [@HochbruckOstermann10] schemes are implemented for a selection of models.

In the same spirit, while the name of the package clearly indicates a disposition towards unidimensional waves, this is by no means an absolute restriction of the package architecture and the two-dimensional shallow-water system has been implemented.

# Features

The package provides users with built-in analysis tools, such as the evaluation of conserved quantities (mass, momentum, energy) for monitoring the accuracy of the discrete approximation.

Visualization tools built on plot recipes of the standard package `Plots.jl` provide easy access to illustrations and preliminary investigation. In particular, users have immediate access to the computed solutions at given time either at computed collocation points or interpolated on any specified grids, as well as coefficients in the frequency domain.

The package implements offers convenient tools to construct initial data; for instance analytic formula for (periodic or solitary) traveling waves of the Serre–Green–Naghdi equations and numerically constructed solutions of the Whitham, Whitham–Boussinesq and Whitham–Green–Naghdi equations, or random periodic initial with prescribed regularity and wavelength.

The package offers a tool for saving (raw) data to a local file, which can be loaded later for future analyses.


# Some implemented models

We name some key models implemented in `WaterWaves1D.jl`, while we refer to the package documentation [@DucheneNavaro] for a complete list together with explicit display of the equations. 

## The water-waves equations

One of the key models is the **water-waves** equations (that is the incompressible Euler equations with a free surface, assuming irrotational flows).
This is the "master equations" from which all other models are derived. 
We use the approach of @DyachenkoKuznetsovSpectorEtAl96 that allows to rewrite the system using only nonlinear pointwise operations or Fourier multipliers. Documented drawbacks of the strategy (see *e.g.* @AmbroseCamassaMarzuolaEtAl22) are that the physical quantities are by design computed on a non-uniform grid of collocation points, and that the method suffers from  anti-resolution for steep waves. Concerning the first issue, we implement interpolation tools for comparison with other models that are solved on uniform grids.

## Shallow-water models

Shallow-water models are derived from the assumption of a small layer depth with respect to the typical horizontal wavelength of the flow. In addition to the standard **shallow-water** system [@Saint-Venant1871], the package proposes among others a family of **Boussinesq** systems [@BonaChenSaut02] as well as fully dispersive **Boussinesq–Whitham** models [@DinvayDutykhKalisch19], the **Serre–Green–Naghdi equations** [@Serre53; @GreenNaghdi76] (and a fully dispersive analogue), and the **Isobe–Kakinuma** model [@Isobe94; @Kakinuma01].

The package provides also provides an educated implementation of a strategy for efficiently solving the Serre–Green–Naghdi equations based on the reformulation of the equations as hyperbolic equations with a constraint and suitable approximation through constraint-relaxation, as proposed by several authors [@FavrieGavrilyuk17; @EscalanteDumbserCastro19; @GuermondPopovTovarEtAl19].


## Small-steepness models

Small-steepness models buikd on the assumption that the amplitude of the wave is small with respect to its horizontal wavelength. They also apply to deep-water or even infinite-layer configurations. Small-steepness models include the **Matsuno** [@Matsuno92] and **Akers-Nicholls** [@AkersNicholls10] systems, as well as the **high-order spectral** model hierarchy (implemented up to quartic nonlinearity) proposed in [@DommermuthYue87; @WestBruecknerJandaEtAl87; @CraigSulem93].
It has been observed that high-order spectral systems suffers from spurious high-frequency amplification, and the code implements a regularization proposed in @DucheneMelinand24.


## Unidirectional models

Standard unidirectional models such as the Korteweg–de Vries, Benjamin–Bona–Mahony and Whitham equations are typically derived assuming unidirectional propagation. They can however be used as building blocks for the propagation of general (unidimensional) waves through the superposition of two counterpropagating waves, assuming that initial data are sufficiently localized [@Lannes; @Emerald21]. This is the 



# State of the field

A multitude of softwares for water wave or ocean dynamics have been developed over many years, and we cannot aim at exhaustiveness in the following description. We focus below on recent (a few years old) related Julia packages, emphasizing scope differences. The main specificity of `WaterWaves1D.jl` with respect to existing alternatives is the desire to offer a common framework for different models in view of facilitating their comparison.

In that respect the most closely related package is arguably  `DispersiveShallowWater.jl` [@LampertWittensteinRanocha25]
which implements structure-preserving numerical methods for dispersive shallow water models, including  unidirectional models such as the Korteweg–de Vries and the Benjamin–Bona–Mahony equations, as well as Boussinesq systems, the Serre–Green–Naghdi equations and its augmented "hyperbolized" system. Note that the implemented numerical methods include the Fourier spectral method but that for the purpose of structure preserving, Fast Fourier transform is avoided and as a consequence the Fourier spectral method is less competitive compared with other methods involving sparse matrices. 


Building on the numerical simulation framework for general conservation laws `Trixi.jl` [@Schlottke-LakemperGassnerRanochaEtAl25], `TrixiShallowWater.jl` [@WintersErsingRanochaEtAl25] provides a discontinuous Galerkin discretization of shallow water equations, including wetting and drying capabilities, sediment transport and variable densities through multilayer equations. While the standard shallow water equations with flat bottom is implemented in `WaterWaves1D.jl`, its primary interest lies in dispersive equations.

`SpectralWaves.jl` [@Paprota25] implements a model for the propagation of water waves over topography that has the same structure as the Isobe–Kakinuma model yet constructed through an expansion akin to high-order spectral models, using the Fourier spectral method.

`GeophysicalFlows.jl` [@ConstantinouWagnerSiegelmanEtAl21] is a package for geophysical fluid dynamics (2D periodic incompressible Navier-Stokes equation and quasi-geostrophic equations) that leverages the pseudo-spectral solver modules of the package `FourierFlows.jl` [@ConstantinouWagnerPaloczyEtAl25]. In Waterwaves1D.jl we have avoided using modules of `FourierFlows.jl` and rather built analogous ones in order not to restrict exclusively to pseudo-spectral methods.

We finally mention the package  `Oceananigans.jl` [@RamadhanWagnerHillEtAl20] that implements finite volume algorithms optimized for high-resolution simulations on GPUs for geophysical models adapted to large-scale ocean dynamics.

# Research impact statement

Features of `WaterWaves1D.jl` were explicitly used for different needs in the following works.

- In @MM4WW the different aspects of the package and of the models are discussed, and (static or dynamic) illustrations are provided using the package;
- In @DucheneKlein22 the numerical strategy is used to investigate numerically open questions such as the stability of traveling waves, and possible finite-time singularities;
- @DucheneMelinand24 uses the package to investigate a mechanism of high-frequency amplification and a proposed numerical remedy, validating the theoretical results obtained therein;
- In @DucheneMarstrander26 the Fourier spectral method applied to hyperbolic systems including the shallow-water system is studied. The package allows to validate some of the theoretical results, and to investigate open problems.
- In @DebusscheMeminMoneyron25, Debussche, Mémin, and Moneyron leverage some of the models of the package to a stochastic framework.


# AI usage disclosure

No generative AI tools were used in the development of this software and its documentation, or the writing of this manuscript, or the preparation of supporting materials.



# References
