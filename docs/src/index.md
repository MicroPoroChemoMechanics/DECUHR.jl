# DECUHR.jl

Pure-Julia port of the **DECUHR** algorithm [EspelidGenz1994](@cite) for automatic
adaptive integration of functions with **homogeneous vertex singularities** over
hyper-rectangular regions.

The package exposes a single algorithm type, [`DecuhrAlgorithm`](@ref), that plugs
directly into the [Integrals.jl](https://github.com/SciML/Integrals.jl) / SciML
ecosystem through the standard `SciMLBase.AbstractIntegralAlgorithm` interface.

## Installation

```julia
import Pkg; Pkg.add("DECUHR")
```

## Quick start

```julia
using Integrals, DECUHR

# ∫₀¹∫₀¹ (x·y)^{-1/2} dx dy = 4
f    = (u, p) -> (u[1] * u[2])^(-0.5)
prob = IntegralProblem(f, (zeros(2), ones(2)))
sol  = solve(prob, DecuhrAlgorithm(singul=2, alpha=-0.5); abstol=1e-8)

println(sol.u)       # ≈ 4.0
println(sol.retcode) # Success
```

## Features

- **Vertex singularities** of the form ``f(\mathbf{x}) \sim \prod_{i=1}^{s} x_i^{\alpha}``
  at the lower-left corner of the domain.
- **Automatic estimation** of the singularity exponent ``\alpha`` and detection of
  logarithmic factors (when `alpha ≤ -singul`, the default).
- **Vector-valued integrands** (`NUMFUN > 1`): all components are integrated
  simultaneously.
- **Four built-in integration rules**: degree-13 (2D, 65 pts), degree-11 (3D,
  127 pts), degree-9 (nD), degree-7 (nD); auto-selected by default.
- **Generic value type** `TV`: works with dual numbers (ForwardDiff) for
  automatic differentiation with respect to parameters `p`.

The algorithm, its rules and its extrapolation are described on the
[Algorithm](algorithm.md) page; the works it rests on are listed on the
[References](@ref references) page.

## License and upstream notice

The Julia port is MIT-licensed. Because it is a translation and modification
of the original Fortran 77 DECUHR routines, the upstream copyright notice of
Espelid & Genz **must be preserved in every copy and every derivative work**.
See the [License](license.md) page for the full text of both the MIT License
and the upstream notice.

DECUHR.jl is an **independent** port: it is not affiliated with, endorsed by, or
supported by the authors of the original DECUHR routines. Please report problems
with this package here rather than to them. The package logo is an original work
made for this port.
