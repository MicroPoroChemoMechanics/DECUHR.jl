# Changelog

## v0.1.3 — 2026-06-01

### Changed

- **Default parameter alignment with Fortran reference:** the default value of
  `wrksub` (maximum number of stored subregions) is increased from 5000 to
  50000 to match the original DECUHR Fortran implementation. This provides
  sufficient refinement budget by default for challenging 3-D singularities and
  other difficult cases. Users requiring smaller memory footprints may still
  pass `wrksub` explicitly.

### Documentation

- Updated examples and API documentation to reflect `wrksub=50000` as the new
  default.

## v0.1.2 — 2026-06-01

### Added

- `sol.stats` now exposes `numevals` (integrand-evaluation count) and `ifail`
  (raw DECUHR return code), useful in particular to interpret a `MaxIters`
  return code.
- Automatic differentiation through `solve` now works **even when `alpha` is
  auto-estimated**: the exponent is estimated on the primal integrand (it is a
  structural property of the singularity, independent of the differentiation
  seed) and the integration then proceeds in the dual-number type. Previously
  this combination returned `ifail = 14`.

### Changed

- **Performance:** `_fully_symmetric_sum!` (Fortran `DEFSHR`) no longer
  allocates a fresh generator copy and evaluation-point vector on every call;
  two caller-preallocated scratch buffers are threaded through instead. The
  result is **bit-for-bit identical**; allocations on a heavy singular case
  (`∫(x·y)^(-1/2)`) drop from ~4.4 MB to ~1.5 MB (−67 %).

### Documentation

- Documented DECUHR's deliberately conservative error estimator: a `MaxIters`
  return code often accompanies a result that is already accurate to better
  than the requested tolerance.

### Validation

- The port was validated **bit-for-bit against the original Fortran `DECUHR`**
  (Espelid & Genz, 1994), compiled locally with gfortran, over six test cases
  with identical parameters: 5 of 6 results match to all 16 significant digits
  (the sixth differs by ~1e-10, pure floating-point reassociation), and the
  `IFAIL` / evaluation counts match on all six. In particular the conservative
  error estimator (a `MaxIters` return code on an already-accurate result) is
  confirmed to be **intrinsic to the algorithm** — the Fortran returns the same
  `IFAIL` — and is therefore documented rather than altered.

### Internal

- Clarified the `BETA`-index clamp in the extrapolation weighted-correction
  branch (a latent off-by-one out-of-bounds in the original Fortran at the
  boundary `UPDATE = N-EMAX`) and hoisted the loop-invariant index.
- Test suite expanded from 9 to 27 assertions (radial singularity, polynomial
  exactness, rule-key sweep, parameter-validation return codes, solution stats,
  ForwardDiff with explicit and auto-estimated `alpha`, small-`emax`
  weighted-correction branch).

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
