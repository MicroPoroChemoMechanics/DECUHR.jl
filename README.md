<p>
  <img src="./docs/src/assets/logo.svg" width="100">
</p>

# DECUHR

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://MicroPoroChemoMechanics.github.io/DECUHR.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://MicroPoroChemoMechanics.github.io/DECUHR.jl/dev/)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/MicroPoroChemoMechanics/DECUHR.jl/blob/main/LICENSE)
[![code style: runic](https://img.shields.io/badge/code_style-%E1%9A%B1%E1%9A%A2%E1%9A%BE%E1%9B%81%E1%9A%B2-pink)](https://github.com/fredrikekre/Runic.jl)

`DECUHR.jl` is a pure-Julia port of the DECUHR algorithm (Espelid &
Genz, 1994) for automatic adaptive integration of functions with
**vertex singularities** over hyper-rectangular regions. It is exposed
as a pluggable algorithm for the
[Integrals.jl](https://docs.sciml.ai/Integrals/stable/) solver stack
via the `SciMLBase.AbstractIntegralAlgorithm` interface.

## Features

- 2-D and 3-D integration on hyper-rectangles.
- Vertex singularity handling with user-supplied or auto-estimated
  exponent ``\alpha`` (singular strength).
- Logarithmic singularities (`logf = k` for ``(\log)^k`` weights).
- Vector-valued integrands (any `NUMFUN`).
- Richardson extrapolation on sub-region averages.
- Reports a `retcode` compatible with Integrals.jl (`Success`,
  `MaxIters`, …).

## Installation

`DECUHR.jl` is released through the dedicated
[MPCM-Registry](https://github.com/MicroPoroChemoMechanics/MPCM-Registry).
Add the registry once, then install the package:

```julia
julia> using Pkg
pkg> registry add https://github.com/MicroPoroChemoMechanics/MPCM-Registry
pkg> add DECUHR
```

## Quick start

```julia
using Integrals, DECUHR

# ∫₀¹∫₀¹ (x·y)^(-0.5) dx dy = 4
f = (u, _) -> (u[1] * u[2])^(-0.5)
prob = IntegralProblem(f, (zeros(2), ones(2)))
sol  = solve(prob, DecuhrAlgorithm(singul = 2, alpha = -0.5); abstol = 1e-8)

@show sol.u           # ≈ 4.0
@show sol.retcode     # Success
```

A broader set of examples is in
[`examples/basic_usage.jl`](./examples/basic_usage.jl):

```shell
julia --project=. examples/basic_usage.jl
```

## Tests

```shell
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Documentation

```shell
julia --project=docs -e 'using Pkg; Pkg.develop(path = "."); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

HTML output is placed in `docs/build/`.

## References

- T.O. Espelid and A. Genz, *DECUHR: An Algorithm for Automatic
  Integration of Singular Functions over a Hyperrectangular Region*,
  Numerical Algorithms **8** (1994) 201–220.
- T.O. Espelid, *On integrating Vertex Singularities using
  Extrapolation*, BIT **34** (1994) 62–79.

## License

The Julia port (this package) is released under the MIT License — see
[LICENSE](LICENSE) for details.

### Upstream / third-party notice

`DECUHR.jl` is a translation and modification of the Fortran 77 DECUHR
routines of Espelid & Genz (Numerical Algorithms 8, 1994). The upstream
distribution carries its own copyright notice, which **must be preserved in
every copy and every derivative work** of this package. That notice is
reproduced verbatim in [NOTICE](NOTICE); redistributors MUST ship both
[LICENSE](LICENSE) and [NOTICE](NOTICE) unmodified.

## Citation

See [CITATION.cff](CITATION.cff) for citation details.

**BibTeX entry:**

```bibtex
@software{decuhr_jl,
  author = {Barthélémy, Jean-François},
  title  = {DECUHR.jl: Adaptive cubature for vertex singularities},
  url    = {https://github.com/MicroPoroChemoMechanics/DECUHR.jl},
  year   = {2026}
}
```

## Credits and Acknowledgements

Developed by [Jean-François Barthélémy](https://github.com/jfbarthelemy),
researcher at [Cerema](https://www.cerema.fr/en) in the research team
[UMR MCD](https://mcd.univ-gustave-eiffel.fr/).

The Fortran 77 → Julia translation and parts of the subsequent codebase
were developed with the assistance of Anthropic's *Claude Code*, under
the author's review and numerical validation.
