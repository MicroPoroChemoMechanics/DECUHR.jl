```@meta
EditURL = "../../../examples/basic_usage.jl"
```

# Examples

The DECUHR algorithm as a pluggable [Integrals.jl](https://github.com/SciML/Integrals.jl)
solver, on integrands whose exact value is known — vertex singularities in two
and three dimensions, a logarithmic singularity, vector-valued and
parametrized integrands, forward-mode differentiation through the quadrature,
and what the return codes and diagnostics actually mean.

This file is the single source of both the documentation page you may be
reading and a standalone script: run it with

    julia --project=. examples/basic_usage.jl

or, from a REPL, `include("examples/basic_usage.jl")`.

````@example examples
using Integrals
using DECUHR
using Printf
````

## 1 — Vertex singularity 2D, known ``\alpha``

Analytical integral:

```math
I = \int_0^1\!\int_0^1 (x_1 x_2)^{-1/2}\, dx_1\, dx_2
  = \left(\int_0^1 x^{-1/2}\,dx\right)^2 = 4
```

````@example examples
f = (u, p) -> (u[1] * u[2])^(-0.5)
prob = IntegralProblem(f, (zeros(2), ones(2)))
sol = solve(prob, DecuhrAlgorithm(singul = 2, alpha = -0.5); abstol = 1.0e-8)

println("I  ≈ ", sol.u)
println("|err| = ", abs(sol.u - 4.0))
println("retcode : ", sol.retcode)
````

## 2 — Automatic estimation of ``\alpha``

Same integral, but without providing ``\alpha``: DECUHR estimates it via DECALP.

````@example examples
sol2 = solve(prob, DecuhrAlgorithm(singul = 2); abstol = 1.0e-7)

println("I  ≈ ", sol2.u)
println("|err| = ", abs(sol2.u - 4.0))
````

## 3 — Smooth integrand (no singularity)

```math
I = \int_0^{\pi/2}\!\int_0^{\pi/2} \sin(x_1)\cos(x_2)\, dx_1\, dx_2 = 1
```

````@example examples
f3 = (u, p) -> sin(u[1]) * cos(u[2])
prob3 = IntegralProblem(f3, (zeros(2), fill(π / 2, 2)))
sol3 = solve(prob3, DecuhrAlgorithm(); abstol = 1.0e-10)

println("I  ≈ ", sol3.u)
println("|err| = ", abs(sol3.u - 1.0))
````

## 4 — Vector-valued integrand (`NUMFUN = 2`)

Two components integrated simultaneously:

```math
I_1 = \int_0^1\!\int_0^1 (x_1^2 + x_2^2)\, dx = \tfrac{2}{3},
\qquad
I_2 = \int_0^1\!\int_0^1 x_1 x_2\, dx = \tfrac{1}{4}
```

````@example examples
f4 = (u, p) -> [u[1]^2 + u[2]^2, u[1] * u[2]]
prob4 = IntegralProblem(f4, (zeros(2), ones(2)))
sol4 = solve(prob4, DecuhrAlgorithm(); abstol = 1.0e-9)

exact4 = [2 / 3, 1 / 4]
println("I  ≈ ", sol4.u)
println("|err| = ", abs.(sol4.u .- exact4))
````

## 5 — 3D singularity

```math
I = \int_0^1\!\int_0^1\!\int_0^1 (x_1 x_2 x_3)^{-1/3}\, dx
  = \left(\frac{3}{2}\right)^3 = \frac{27}{8}
```

The 3-D vertex singularity is challenging; `wrksub=50000` (now the default)
allows sufficient refinement:

````@example examples
f5 = (u, p) -> (u[1] * u[2] * u[3])^(-1 / 3)
prob5 = IntegralProblem(f5, (zeros(3), ones(3)))
sol5 = solve(
    prob5, DecuhrAlgorithm(singul = 3, alpha = -1 / 3);
    abstol = 1.0e-7, reltol = 1.0e-7, maxiters = 1_500_000
)

println("I  ≈ ", sol5.u)
println("|err| = ", abs(sol5.u - (3 / 2)^3))
````

## 6 — Logarithmic singularity

```math
I = \int_0^1\!\int_0^1 -\log(x_1 x_2)\, dx = 2
```

````@example examples
f6 = (u, p) -> -log(u[1] * u[2])
prob6 = IntegralProblem(f6, (zeros(2), ones(2)))
sol6 = solve(prob6, DecuhrAlgorithm(singul = 2, alpha = 0.0, logf = 1); abstol = 1.0e-8)

println("I  ≈ ", sol6.u)
println("|err| = ", abs(sol6.u - 2.0))
````

## 7 — Parametrized integral

Integral depending on a parameter ``\lambda`` passed via `p`:

```math
I(\lambda) = \int_0^1\!\int_0^1
(x_1 x_2)^{-1/2}\, e^{-\lambda(x_1 + x_2)}\, dx
```

For ``\lambda = 0``: ``I(0) = 4``.

````@example examples
f7 = (u, p) -> (u[1] * u[2])^(-0.5) * exp(-p[1] * (u[1] + u[2]))
prob7 = IntegralProblem(f7, (zeros(2), ones(2)), [0.0])

for λ in (0.0, 0.5, 1.0, 2.0)
    s = solve(
        remake(prob7, p = [λ]),
        DecuhrAlgorithm(singul = 2, alpha = -0.5);
        abstol = 1.0e-8
    )
    @printf "λ = %.1f  →  I ≈ %.6f\n" λ s.u
end
````

## 8 — Automatic differentiation with ForwardDiff

The integrand is parametrized by ``\lambda``. We compute ``dI/d\lambda`` and
``d^2I/d\lambda^2`` in forward AD mode, without finite differences.

```math
I(\lambda)
= \int_0^1\!\int_0^1 (x_1 x_2)^{-1/2}\, e^{-\lambda(x_1+x_2)}\, dx
```

**Analytical derivative:**

```math
\frac{dI}{d\lambda}
= -\int_0^1\!\int_0^1 (x_1+x_2)\,(x_1 x_2)^{-1/2}\,
    e^{-\lambda(x_1+x_2)}\, dx
```

At ``\lambda = 0``:

```math
\frac{dI}{d\lambda}\bigg|_{\lambda=0}
= -2\int_0^1 x^{-1/2}\,dx \cdot \int_0^1 x^{1/2}\,dx
= -2 \cdot 2 \cdot \tfrac{2}{3} = -\tfrac{8}{3}
```

````@example examples
using ForwardDiff

f8 = (u, p) -> (u[1] * u[2])^(-0.5) * exp(-p[1] * (u[1] + u[2]))
prob8 = IntegralProblem(f8, (zeros(2), ones(2)), [0.0])

# Function I(λ). Here λ does not change the singularity structure, so we may
# either supply alpha explicitly or let it be auto-estimated — both differentiate.
I(λ) = solve(
    remake(prob8, p = [λ]),
    DecuhrAlgorithm(singul = 2, alpha = -0.5);
    abstol = 1.0e-7
).u

# First derivative
dI = ForwardDiff.derivative(I, 0.0)

# Second derivative (second-order nesting)
d2I = ForwardDiff.derivative(λ -> ForwardDiff.derivative(I, λ), 0.0)

exact_dI = -8 / 3
println("dI/dλ  ≈ ", dI, "  (exact = ", exact_dI, ")")
println("|err|    = ", abs(dI - exact_dI))
println("d²I/dλ² ≈ ", d2I)
````

The same derivative is obtained **without** supplying `alpha`: it is
auto-estimated on the primal integrand (a structural property of the
singularity, independent of the differentiation seed) and the integration then
runs in the dual-number type.

````@example examples
Iauto(λ) = solve(
    remake(prob8, p = [λ]),
    DecuhrAlgorithm(singul = 2);   # alpha auto-estimated
    abstol = 1.0e-7, maxiters = 300_000
).u

println(
    "dI/dλ (auto-α) ≈ ", ForwardDiff.derivative(Iauto, 0.0),
    "   (exact = ", -8 / 3, ")"
)
````

### Multi-parameter gradient

We add a second parameter ``\mu`` controlling the singularity exponent:

```math
I(\lambda, \mu)
= \int_0^1\!\int_0^1 (x_1 x_2)^{\mu}\, e^{-\lambda(x_1+x_2)}\, dx,
\qquad \mu > -1
```

**Analytical values at** ``(\lambda, \mu) = (0,\,-0.3)``:

```math
I = \frac{1}{(1+\mu)^2},
\quad
\frac{\partial I}{\partial\lambda} = -\frac{2}{(1+\mu)(2+\mu)},
\quad
\frac{\partial I}{\partial\mu} = \frac{-2}{(1+\mu)^3}
```

````@example examples
f9 = (u, p) -> (u[1] * u[2])^p[2] * exp(-p[1] * (u[1] + u[2]))
prob9 = IntegralProblem(f9, (zeros(2), ones(2)), [0.0, -0.3])

# I as a function of p = [λ, μ] — singularity exponent = p[2], held fixed
Ivec(p) = solve(
    remake(prob9, p = p),
    DecuhrAlgorithm(singul = 2, alpha = -0.3);   # alpha fixed at μ₀
    abstol = 1.0e-7
).u

p0 = [0.0, -0.3]
grad = ForwardDiff.gradient(Ivec, p0)

μ = p0[2]
exact_I = 1 / (1 + μ)^2
# dI/dλ|λ=0 = -∫∫ (x₁+x₂)(x₁x₂)^μ dx = -2/((μ+1)(μ+2))
exact_dIdλ = -2 / ((1 + μ) * (2 + μ))
# dI/dμ = d/dμ [1/(1+μ)²] = -2/(1+μ)³
exact_dIdμ = -2 / (1 + μ)^3

println("I  ≈ ", exact_I)
println("∇I ≈ ", grad)
println("exact ∇I = [", exact_dIdλ, ", ", exact_dIdμ, "]")
println("|err|    = ", abs.(grad .- [exact_dIdλ, exact_dIdμ]))
````

!!! warning "Alpha held fixed during gradient computation"
    When differentiating with respect to ``\mu`` (the singularity exponent),
    `alpha` in `DecuhrAlgorithm` must remain a constant `Float64` —
    it controls the extrapolation rule, not the value of the integral.
    For an exact gradient in ``\mu``, one must evaluate at ``\mu_0``
    and accept that the quadrature error introduces a bias of order
    ``O(\text{abstol})``.

## 9 — Budget control

What happens when the budget cannot meet the tolerance depends on *how* short
it is, and the two cases are not the same. With a budget too small to complete
even the first subdivision of a 3-D vertex singularity, the solve fails
outright and there is **no** estimate to salvage — `u` comes back as zero, not
as a poor approximation:

````@example examples
sol9 = solve(
    prob5,
    DecuhrAlgorithm(singul = 3, alpha = -1 / 3);
    abstol = 1.0e-14,   # very tight tolerance
    maxiters = 500      # far too small a budget
)

println("retcode : ", sol9.retcode)
println("u       : ", sol9.u)
println("resid   : ", sol9.resid)
println("ifail   : ", sol9.stats.ifail)
````

Always test the return code before using `u`. With a budget large enough to
refine, but still short of the requested tolerance, the return code becomes
`MaxIters` and the value *is* usable — accurate to about 1 %, which is not the
requested `1e-14` but is a genuine estimate:

````@example examples
sol9b = solve(
    prob5,
    DecuhrAlgorithm(singul = 3, alpha = -1 / 3);
    abstol = 1.0e-14,
    maxiters = 50_000
)

println("retcode : ", sol9b.retcode)
println("best estimate : ", sol9b.u)
println("|err| ≈ ", abs(sol9b.u - (3 / 2)^3))
````

## 10 — Diagnostics: `sol.stats` and `sol.resid`

`sol.stats` reports the number of integrand evaluations and the raw DECUHR
code, and `sol.resid` carries the algorithm's own error estimate. Asking the
2-D vertex singularity for `1e-12` exhausts the budget: the return code is
`MaxIters` and the value is accurate to about `1e-4`.

Note what `resid` says next to the true error below. The estimate is not a
guaranteed bound: here it is smaller than the actual error by an order of
magnitude. Treat it as an indication of how far the extrapolation has
converged, not as a certificate — and where the value matters, confirm it by
tightening `maxiters` until `u` stops moving.

````@example examples
sol10 = solve(prob, DecuhrAlgorithm(singul = 2, alpha = -0.5); abstol = 1.0e-12)

println("retcode   : ", sol10.retcode)
println("numevals  : ", sol10.stats.numevals) # integrand evaluations
println("ifail     : ", sol10.stats.ifail)    # raw DECUHR code
println("resid     : ", sol10.resid)          # the algorithm's own estimate
println("I ≈ ", sol10.u, "   true |err| = ", abs(sol10.u - 4.0))
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

