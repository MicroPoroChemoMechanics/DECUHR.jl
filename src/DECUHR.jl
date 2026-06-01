# ***************************************************************************
# * All the software  contained in this library  is protected by copyright. *
# * Permission  to use, copy, modify, and  distribute this software for any *
# * purpose without fee is hereby granted, provided that this entire notice *
# * is included  in all copies  of any software which is or includes a copy *
# * or modification  of this software  and in all copies  of the supporting *
# * documentation for such software.                                        *
# ***************************************************************************
# * THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED *
# * WARRANTY. IN NO EVENT, NEITHER  THE AUTHORS, NOR THE PUBLISHER, NOR ANY *
# * MEMBER  OF THE EDITORIAL BOARD OF  THE JOURNAL  "NUMERICAL ALGORITHMS", *
# * NOR ITS EDITOR-IN-CHIEF, BE  LIABLE FOR ANY ERROR  IN THE SOFTWARE, ANY *
# * MISUSE  OF IT  OR ANY DAMAGE ARISING OUT OF ITS USE. THE ENTIRE RISK OF *
# * USING THE SOFTWARE LIES WITH THE PARTY DOING SO.                        *
# ***************************************************************************
# * ANY USE  OF THE SOFTWARE  CONSTITUTES  ACCEPTANCE  OF THE TERMS  OF THE *
# * ABOVE STATEMENT.                                                        *
# ***************************************************************************
#
# Reference: T.O. Espelid and A. Genz, "DECUHR: An Algorithm for Automatic
# Integration of Singular Functions over a Hyperrectangular Region",
# Numerical Algorithms 8 (1994), pp. 201-220.
#
"""
    DECUHR

Pure-Julia port of the DECUHR algorithm (Espelid & Genz, 1994) for automatic
integration of functions with vertex singularities over hyper-rectangular
regions.

Exposed as a pluggable algorithm for the SciML **Integrals.jl** ecosystem
through the `SciMLBase.AbstractIntegralAlgorithm` interface.

## Quick start

```julia
using Integrals, DECUHR

# ∫∫₀¹ (x·y)^{-0.5} dx dy = π²/4
f = (u, p) -> (u[1]*u[2])^(-0.5)
prob = IntegralProblem(f, (zeros(2), ones(2)))
sol  = solve(prob, DecuhrAlgorithm(singul=2, alpha=-0.5); abstol=1e-8)
```

## References

- T.O. Espelid and A. Genz, *DECUHR: An Algorithm for Automatic Integration
  of Singular Functions over a Hyperrectangular Region*,
  Numerical Algorithms 8 (1994), pp. 201-220.
- T.O. Espelid, *On integrating Vertex Singularities using Extrapolation*,
  BIT 34 (1994), pp. 62-79.

## License

This package is released under the MIT License (`LICENSE`). As a port and
modification of the Fortran 77 DECUHR routines, it additionally carries the
upstream copyright notice of Espelid & Genz, reproduced verbatim in `NOTICE`.
Both files MUST be shipped with every copy and every derivative work.
"""
module DECUHR

using SciMLBase, Integrals

include("rules.jl")
include("extrapolation.jl")
include("alpha_estimation.jl")
include("adaptive.jl")

export DecuhrAlgorithm

# ============================================================
# Algorithm struct
# ============================================================
"""
    DecuhrAlgorithm(; key=0, singul=1, alpha=-2.0, logf=0,
                      wrksub=50000, emax=20, minpts=0)

DECUHR adaptive integration algorithm for functions with vertex singularities.

## Keyword arguments

| Argument | Default | Description |
|:---------|:--------|:------------|
| `key`    | `0`     | Rule selector: 0=auto, 1=2D deg-13 (65 pts), 2=3D deg-11 (127 pts), 3=deg-9 (nD), 4=deg-7 (nD) |
| `singul` | `1`     | Dimension of the singularity: the vertex `(a[1],…,a[singul])` is singular |
| `alpha`  | `-2.0`  | Exponent of the homogeneous singularity. Auto-estimated when `alpha ≤ -singul` |
| `logf`   | `0`     | `1` if a logarithmic factor is present; `0` otherwise (auto-detected when `alpha ≤ -singul`) |
| `wrksub` | `50000` | Maximum number of stored subregions |
| `emax`   | `20`    | Maximum number of Richardson extrapolation steps |
| `minpts` | `0`     | Minimum number of integrand evaluations |

## Notes

- The singularity must be located at the **lower-left vertex** `a[1:singul]`
  of the integration domain.  Shift the domain if necessary.
- `alpha` must satisfy `alpha > -singul`.  Setting `alpha ≤ -singul` (the
  default) triggers automatic estimation via `DECALP`.
- For regular (non-singular) integrands set `singul=1` and leave `alpha` at
  its default; the algorithm degrades gracefully to ordinary adaptive Gauss.

## Return codes

| `sol.retcode`         | Meaning |
|:----------------------|:--------|
| `ReturnCode.Success`  | All error estimates meet the requested tolerances |
| `ReturnCode.MaxIters` | Budget (`maxiters`) exhausted before convergence |
| `ReturnCode.Failure`  | Invalid input parameters (see `sol.stats`) |

## Interpreting the result and `MaxIters`

DECUHR's error estimator is deliberately **conservative** (it inflates the
pure extrapolation error by a heuristic factor, exactly as in the original
Fortran `DEXTHR`).  As a consequence a returned `retcode = MaxIters` does **not**
mean the result is wrong: the value in `sol.u` is frequently accurate to far
better than the requested tolerance even when the *estimated* error
(`sol.resid`) has not dropped below it.  When this happens, either raise
`maxiters`/`wrksub` to let the estimate catch up, or simply trust `sol.u` and
inspect `sol.resid`.

`sol.stats` exposes diagnostics:

| Field                | Meaning |
|:---------------------|:--------|
| `sol.stats.numevals` | Number of integrand evaluations performed |
| `sol.stats.ifail`    | Raw DECUHR `IFAIL` code (0 = success, 1 = budget hit, …) |

## Automatic differentiation

`solve` is differentiable through ForwardDiff with respect to integrand
parameters, **including** when `alpha` is auto-estimated: the exponent is
estimated on the primal integrand (it is a structural property of the
singularity) and the integration is then carried out in the dual number type.
"""
struct DecuhrAlgorithm <: SciMLBase.AbstractIntegralAlgorithm
    key::Int
    singul::Int
    alpha::Float64
    logf::Int
    wrksub::Int
    emax::Int
    minpts::Int
end

function DecuhrAlgorithm(;
        key = 0,
        singul = 1,
        alpha = -2.0,   # ≤ -singul triggers auto-estimation
        logf = 0,
        wrksub = 50000,
        emax = 20,
        minpts = 0
    )
    return DecuhrAlgorithm(key, singul, alpha, logf, wrksub, emax, minpts)
end

# ============================================================
# Internal driver — mirrors Fortran DECUHR
# ============================================================
const _MAXDIM = 15

"""
    _primal_to_float(x) -> Float64

Return the primal `Float64` value of a real number, peeling off any
ForwardDiff-style dual layers.  DECUHR does not depend on ForwardDiff, so we
duck-type: a dual number exposes a `value` field holding its primal part
(possibly itself a dual), which we recurse into.  Plain reals convert directly.

This lets the `α`/`logf` auto-estimation (`DECALP`, Float64-only) run on the
*primal* integrand when the user differentiates through `solve` — `α` is the
structural degree of the vertex singularity and does not itself depend on the
differentiation seed, so estimating it on the primal values is exact.
"""
_primal_to_float(x::AbstractFloat) = Float64(x)
function _primal_to_float(x::Real)
    return hasfield(typeof(x), :value) ? _primal_to_float(getfield(x, :value)) : Float64(x)
end

"""
    _decuhr_driver(ndim, numfun, a, b, minpts, maxpts, funsub,
                   singul, alpha_in, logf_in, epsabs, epsrel,
                   key, wrksub, emax)

Low-level driver.  Validates parameters, optionally estimates `alpha`, then
calls `_adaptive_integrate!`.  Returns `(result, abserr, neval, ifail)`.
"""
function _decuhr_driver(
        ndim::Int, numfun::Int,
        a::AbstractVector, b::AbstractVector,
        minpts::Int, maxpts::Int,
        funsub, ::Type{TV},
        singul::Int, alpha_in::Float64, logf_in::Int,
        epsabs::Real, epsrel::Real,
        key::Int, wrksub::Int, emax::Int
    ) where {TV}

    # --- Parameter validation ---
    num, maxsub, keyf, wtleng, ifail = _check_params(
        _MAXDIM, ndim, numfun, a, b, singul, logf_in,
        minpts, maxpts, epsabs, epsrel, key, emax, wrksub
    )

    if ifail != 0
        return zeros(TV, numfun), zeros(TV, numfun), 0, ifail
    end

    alpha = alpha_in
    logf = logf_in

    # --- Auto-estimate ALPHA (DECALP) when alpha ≤ -singul ---
    if alpha <= -Float64(singul)
        alpha_ref = Ref(alpha)
        logf_ref = Ref(logf)
        x_work = zeros(Float64, ndim)
        funs_work = zeros(Float64, numfun)

        # DECALP works in Float64.  For Float64 integrands use `funsub` directly;
        # for non-Float64 value types (e.g. ForwardDiff duals when differentiating
        # through `solve`) estimate α on the *primal* values via `_primal_to_float`.
        # α is the structural homogeneity degree of the singularity and does not
        # depend on the differentiation seed, so this is exact and keeps `solve`
        # differentiable even when α is auto-estimated.
        alpha_funsub = if TV === Float64
            funsub
        else
            funs_tv = zeros(TV, numfun)
            function (x, out64)
                funsub(x, funs_tv)
                @inbounds for j in eachindex(out64)
                    out64[j] = _primal_to_float(funs_tv[j])
                end
                return nothing
            end
        end

        _estimate_alpha!(
            ndim, a, b, alpha_ref, singul, logf_ref,
            alpha_funsub, x_work, funs_work
        )
        alpha = alpha_ref[]
        logf = logf_ref[]

        if alpha <= -Float64(singul)
            return zeros(TV, numfun), zeros(TV, numfun), 0, 14
        end
    end

    # --- Adaptive integration ---
    result, abserr, neval, ifail = _adaptive_integrate!(
        ndim, numfun, a, b, maxsub, funsub, TV,
        singul, alpha, logf,
        Float64(epsabs), Float64(epsrel), keyf, num, wtleng,
        emax, minpts, maxpts
    )

    return result, abserr, neval, ifail
end

# ============================================================
# Integrals.__solvebp_call — Integrals.jl integration hook
# ============================================================
function Integrals.__solvebp_call(
        cache::Integrals.IntegralCache,
        alg::DecuhrAlgorithm,
        sensealg, domain, p;
        abstol = 1.0e-8,
        reltol = 1.0e-6,
        maxiters = 100_000,
        kwargs...
    )

    lb, ub = domain
    ndim = length(lb)
    f = cache.f   # IntegralFunction wrapping the user's integrand

    # Detect numfun and value type TV via a single test evaluation at the midpoint.
    # TV may be a dual-number type when differentiating w.r.t. parameters p.
    xmid = (lb .+ ub) ./ 2
    test_out = f(xmid, p)
    numfun = test_out isa Number ? 1 : length(test_out)
    TV = test_out isa Number ? typeof(test_out) : eltype(test_out)

    # Wrap SciML f(u, p) → DECUHR funsub(x, funvls).
    # funvls may be a SubArray (column-view from rules.jl); no type restriction.
    if numfun == 1
        funsub = (x, funvls) -> (funvls[1] = f(x, p); nothing)
    else
        funsub = (x, funvls) -> (copyto!(funvls, f(x, p)); nothing)
    end

    result, abserr, neval, ifail = _decuhr_driver(
        ndim, numfun,
        lb, ub,
        alg.minpts, maxiters,
        funsub, TV,
        alg.singul, Float64(alg.alpha), alg.logf,
        abstol, reltol,
        alg.key, alg.wrksub, alg.emax
    )

    retcode = if ifail == 0
        SciMLBase.ReturnCode.Success
    elseif ifail == 1
        SciMLBase.ReturnCode.MaxIters
    else
        SciMLBase.ReturnCode.Failure
    end

    u = numfun == 1 ? result[1] : result
    err = numfun == 1 ? abserr[1] : abserr

    prob = Integrals.build_problem(cache)
    # Expose the integrand-evaluation count (and raw IFAIL) via `sol.stats`.
    # Useful in particular to diagnose `MaxIters`: the result may already meet
    # the tolerance even though DECUHR's (conservative) error estimate did not.
    return SciMLBase.build_solution(
        prob, alg, u, err;
        retcode = retcode,
        stats = (numevals = neval, ifail = ifail)
    )
end

end  # module DECUHR
