# Changelog

## v0.1.1 — Integrals.jl tuple-domain API

### Changed

- Documentation, examples and tests now construct `IntegralProblem`
  with the **tuple-domain** form `IntegralProblem(f, (lb, ub))`
  (and `IntegralProblem(f, (lb, ub), p)` for parameterised integrands).
  The legacy three-positional-argument form `IntegralProblem(f, lb, ub)`
  is no longer supported by Integrals.jl / SciMLBase (≥ 5.4): the upper
  bound `ub` was silently swallowed as the parameter `p` and `lb` alone
  became the domain, so `DecuhrAlgorithm` received scalar bounds and
  raised a `BoundsError`. No change to `DecuhrAlgorithm` or the solver
  itself — only the documented call syntax.

### Compat

- `[compat]` lower bound for `Integrals` raised to `"5.4"` (was `"5"`),
  the line that ships the modern tuple-domain / `ChangeOfVariables`
  domain handling this release targets. `SciMLBase` stays `"2, 3"`:
  DECUHR is exercised on both SciMLBase 2.155 (via MeanFieldHom.jl, whose
  `OrdinaryDiffEq 6` pins SciMLBase 2) and SciMLBase 3.16 (standalone).

## v0.1.0 — Initial release

First public release of DECUHR.jl, hosted on Codeberg
(`MicroPoroChemoMechanics/DECUHR.jl`) and registered in MPCM-Registry.

Pure-Julia port of the DECUHR algorithm
([Espelid & Genz, *Numerical Algorithms* 8 (1994), 201–220](https://doi.org/10.1007/BF02142690))
for automatic adaptive integration of functions with **vertex
singularities** over hyper-rectangular regions.

### Features

- `DecuhrAlgorithm` plugs into the SciML
  [Integrals.jl](https://github.com/SciML/Integrals.jl) solver stack
  via `SciMLBase.AbstractIntegralAlgorithm`.
- Genz–Malik cubature rule plus corner-aware variants
  (`src/rules.jl`).
- Richardson extrapolation on sub-region averages
  (`src/extrapolation.jl`).
- Automatic estimation of the singular strength `α` from the empirical
  decay of sub-region contributions (`src/alpha_estimation.jl`); user
  override available via `DecuhrAlgorithm(; alpha = …)`.
- Adaptive subdivision loop (`src/adaptive.jl`) with configurable
  workspace size (`wrksub`), extrapolation order (`emax`) and minimum
  function-evaluation count (`minpts`).
- 2D and 3D integration; `singul ∈ {1, 2, 3}` selects edge / face /
  corner singularity structure.
- Logarithmic singularity factor `logf ∈ {0, 1, …}`
  (e.g. `f(x) ~ x^α · (-log x)^logf`).
- Vector-valued integrand support.
- ForwardDiff-compatible throughout: gradients propagate through both
  integrand and integration domain when the rule is applied to a
  problem parameterised by `Dual` numbers.

### Infrastructure

- Forgejo workflows on Codeberg: CI, Documentation, Release, Runic,
  Zenodo (`.forgejo/workflows/`).
- Multi-version documentation deployment via `docs/deploy_docs.jl`
  (`stable/`, `dev/`, `vX.Y.Z/`).
- Runic.yml is `workflow_dispatch` only to keep auto-format under
  manual control.
- GitHub workflows kept in `.github/workflows/` as a dormant return
  path.

### Provenance

The Julia port is released under the **MIT** license. The upstream
Fortran package by Espelid & Genz (1994) carries its own copyright
notice, reproduced verbatim in `NOTICE` and in
`docs/src/license.md` — both must travel with any redistribution of
this package.
