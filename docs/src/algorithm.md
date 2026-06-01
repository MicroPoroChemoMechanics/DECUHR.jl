# Algorithm

## Problem class

DECUHR targets integrals of the form

```math
I = \int_{\mathbf{a}}^{\mathbf{b}} f(\mathbf{x})\, d\mathbf{x}
```

where ``f`` may have a **homogeneous vertex singularity** at the lower-left
corner ``(a_1, \ldots, a_s)`` of the domain:

```math
f(\mathbf{x}) \;\sim\; g(\mathbf{x})\,
\prod_{i=1}^{s}(x_i - a_i)^{\alpha}
\quad \text{as } \mathbf{x} \to \mathbf{a}_{1:s},
```

where ``\alpha > -s`` and ``g`` is smooth near the corner.
A logarithmic factor ``\log\!\bigl(\prod_i(x_i-a_i)\bigr)`` is also supported.

!!! note "Domain convention"
    The singularity must be located at the **lower-left vertex** of the
    integration domain, i.e. at ``\mathbf{a}_{1:\text{singul}}``.
    Shift the domain beforehand if necessary.

## Strategy

The algorithm combines three ideas:

1. **Adaptive subdivision** (DEADHR / DESBHR).
   The domain is maintained as a pool of subregions stored in a **binary
   max-heap** keyed on the local error estimate.  At each step the
   worst-error region is subdivided:
   - *Singular region* (touching the singularity corner): cut into
     `singul + 1` children by halving each singular dimension in turn.
   - *Regular region*: bisect along the axis with the largest fourth
     difference.

2. **Richardson extrapolation** (DEXTHR).
   Each time the singular region is subdivided, a new term ``U_k`` of a
   geometric series is appended.  The tableau is extrapolated to accelerate
   convergence from ``O(h^\alpha)`` to machine precision.

3. **Automatic singularity estimation** (DECALP).
   When `alpha ≤ -singul` (the default), the exponent is estimated by
   evaluating the integrand along a ray toward the singular vertex and
   extrapolating the log-ratio sequence with Aitken / Bjorstad–Grosse–Dahlquist
   acceleration.

## Integration rules

Four Genz–Malik-style fully-symmetric rules are available, selected via `key`:

| `key` | Rule | Dimensions | Points per region |
|:-----:|------|:----------:|:-----------------:|
| `1`   | Degree-13 | 2D only | 65 |
| `2`   | Degree-11 | 3D only | 127 |
| `3`   | Degree-9  | any nD  | ``1 + 2n + 6n^2 + \ldots + 2^n`` |
| `4`   | Degree-7  | any nD  | ``1 + 2n(n+2) + 2^n`` |
| `0`   | Auto      | —       | key=1 if n=2, key=2 if n=3, else key=3 |

## Generic value type and automatic differentiation

All internal arrays that accumulate integrand values are parameterized on a
value type `TV`, inferred at runtime from a test evaluation:

```julia
TV = typeof(f(midpoint, p))   # e.g. Float64 or ForwardDiff.Dual{...}
```

This means the algorithm works out of the box with **dual numbers** when
differentiating through the integral with respect to parameters `p`:

```julia
using ForwardDiff, Integrals, DECUHR

f = (u, p) -> exp(-p[1] * (u[1]^2 + u[2]^2))
prob = IntegralProblem(f, (zeros(2), ones(2)))

# Derivative of the integral w.r.t. p[1] at p₀ = [1.0]
dI_dp = ForwardDiff.gradient(p -> solve(prob, DecuhrAlgorithm(alpha=0.0);
                                         abstol=1e-8).u, [1.0])
```

!!! warning "Alpha auto-estimation and AD"
    Automatic estimation of ``\alpha`` (triggered when `alpha ≤ -singul`)
    uses `Float64` arithmetic internally and is **not compatible** with
    dual-number value types.  When differentiating, always supply `alpha`
    explicitly:
    ```julia
    DecuhrAlgorithm(singul=2, alpha=-0.5)   # correct
    DecuhrAlgorithm(singul=2)               # will error with dual numbers
    ```

## Error codes (`ifail`)

The low-level driver returns an integer `ifail`:

| Code | Meaning |
|:----:|---------|
| `0`  | Success — all tolerances satisfied |
| `1`  | Budget (`maxiters`) exhausted |
| `2`–`13` | Invalid input parameters (key, ndim, numfun, bounds, budget, tolerances, singul, logf, minpts, emax, wrksub) |
| `14` | Integral may not converge (``\alpha \leq -\text{singul}`` after estimation, or dual-number type without explicit alpha) |
| `15` | `emax < 1` |
| `17` | `wrksub` too small for the given budget |

These are mapped to `SciMLBase.ReturnCode` values by the Integrals.jl integration hook:
`ifail=0` → `Success`, `ifail=1` → `MaxIters`, otherwise → `Failure`.
