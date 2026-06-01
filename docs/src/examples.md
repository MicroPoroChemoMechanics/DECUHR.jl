# Examples

```@setup shared
using Integrals, DECUHR, Printf
```

## 1 — Vertex singularity 2D, known ``\alpha``

Analytical integral:

```math
I = \int_0^1\!\int_0^1 (x_1 x_2)^{-1/2}\, dx_1\, dx_2
  = \left(\int_0^1 x^{-1/2}\,dx\right)^2 = 4
```

```@example shared
f    = (u, p) -> (u[1] * u[2])^(-0.5)
prob = IntegralProblem(f, (zeros(2), ones(2)))
sol  = solve(prob, DecuhrAlgorithm(singul=2, alpha=-0.5); abstol=1e-8)

println("I  ≈ ", sol.u)
println("|err| = ", abs(sol.u - 4.0))
println("retcode : ", sol.retcode)
```

---

## 2 — Automatic estimation of ``\alpha``

Same integral, but without providing ``\alpha``: DECUHR estimates it via DECALP.

```@example shared
sol2 = solve(prob, DecuhrAlgorithm(singul=2); abstol=1e-7)

println("I  ≈ ", sol2.u)
println("|err| = ", abs(sol2.u - 4.0))
```

---

## 3 — Smooth integrand (no singularity)

```math
I = \int_0^{\pi/2}\!\int_0^{\pi/2} \sin(x_1)\cos(x_2)\, dx_1\, dx_2 = 1
```

```@example shared
f3   = (u, p) -> sin(u[1]) * cos(u[2])
prob3 = IntegralProblem(f3, (zeros(2), fill(π/2, 2)))
sol3  = solve(prob3, DecuhrAlgorithm(); abstol=1e-10)

println("I  ≈ ", sol3.u)
println("|err| = ", abs(sol3.u - 1.0))
```

---

## 4 — Vector-valued integrand (`NUMFUN = 2`)

Two components integrated simultaneously:

```math
I_1 = \int_0^1\!\int_0^1 (x_1^2 + x_2^2)\, dx = \tfrac{2}{3},
\qquad
I_2 = \int_0^1\!\int_0^1 x_1 x_2\, dx = \tfrac{1}{4}
```

```@example shared
f4   = (u, p) -> [u[1]^2 + u[2]^2,  u[1]*u[2]]
prob4 = IntegralProblem(f4, (zeros(2), ones(2)))
sol4  = solve(prob4, DecuhrAlgorithm(); abstol=1e-9)

exact4 = [2/3, 1/4]
println("I  ≈ ", sol4.u)
println("|err| = ", abs.(sol4.u .- exact4))
```

---

## 5 — 3D singularity

```math
I = \int_0^1\!\int_0^1\!\int_0^1 (x_1 x_2 x_3)^{-1/3}\, dx
  = \left(\frac{3}{2}\right)^3 = \frac{27}{8}
```

The 3-D vertex singularity is challenging; `wrksub=50000` (now the default)
allows sufficient refinement:

```@example shared
f5   = (u, p) -> (u[1] * u[2] * u[3])^(-1/3)
prob5 = IntegralProblem(f5, (zeros(3), ones(3)))
sol5  = solve(prob5, DecuhrAlgorithm(singul=3, alpha=-1/3);
              abstol=1e-7, reltol=1e-7, maxiters=1_500_000)

println("I  ≈ ", sol5.u)
println("|err| = ", abs(sol5.u - (3/2)^3))
```

---

## 6 — Logarithmic singularity

```math
I = \int_0^1\!\int_0^1 -\log(x_1 x_2)\, dx = 2
```

```@example shared
f6   = (u, p) -> -log(u[1] * u[2])
prob6 = IntegralProblem(f6, (zeros(2), ones(2)))
sol6  = solve(prob6, DecuhrAlgorithm(singul=2, alpha=0.0, logf=1); abstol=1e-8)

println("I  ≈ ", sol6.u)
println("|err| = ", abs(sol6.u - 2.0))
```

---

## 7 — Parametrised integral

Integral depending on a parameter ``\lambda`` passed via `p`:

```math
I(\lambda) = \int_0^1\!\int_0^1
(x_1 x_2)^{-1/2}\, e^{-\lambda(x_1 + x_2)}\, dx
```

For ``\lambda = 0``: ``I(0) = 4``.

```@example shared
f7   = (u, p) -> (u[1] * u[2])^(-0.5) * exp(-p[1] * (u[1] + u[2]))
prob7 = IntegralProblem(f7, (zeros(2), ones(2)), [0.0])

for λ in (0.0, 0.5, 1.0, 2.0)
    s = solve(remake(prob7, p=[λ]),
              DecuhrAlgorithm(singul=2, alpha=-0.5);
              abstol=1e-8)
    @printf "λ = %.1f  →  I ≈ %.6f\n" λ s.u
end
```

---

## 8 — Automatic differentiation with ForwardDiff

The integrand is parametrised by ``\lambda``. We compute
``dI/d\lambda`` and ``d^2I/d\lambda^2`` in forward AD mode,
without finite differences.

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

```@example shared
using ForwardDiff

f8   = (u, p) -> (u[1] * u[2])^(-0.5) * exp(-p[1] * (u[1] + u[2]))
prob8 = IntegralProblem(f8, (zeros(2), ones(2)), [0.0])

# Function I(λ).  Here λ does not change the singularity structure, so we may
# either supply alpha explicitly or let it be auto-estimated — both differentiate.
I(λ) = solve(remake(prob8, p=[λ]),
             DecuhrAlgorithm(singul=2, alpha=-0.5);
             abstol=1e-7).u

# First derivative
dI  = ForwardDiff.derivative(I, 0.0)

# Second derivative (second-order nesting)
d2I = ForwardDiff.derivative(λ -> ForwardDiff.derivative(I, λ), 0.0)

exact_dI  = -8/3
println("dI/dλ  ≈ ", dI,  "  (exact = ", exact_dI, ")")
println("|err|    = ", abs(dI - exact_dI))
println("d²I/dλ² ≈ ", d2I)
```

The same derivative is obtained **without** supplying `alpha`: it is
auto-estimated on the primal integrand (a structural property of the
singularity, independent of the differentiation seed) and the integration then
runs in the dual-number type.

```@example shared
Iauto(λ) = solve(remake(prob8, p=[λ]),
                 DecuhrAlgorithm(singul=2);   # alpha auto-estimated
                 abstol=1e-7, maxiters=300_000).u

println("dI/dλ (auto-α) ≈ ", ForwardDiff.derivative(Iauto, 0.0),
        "   (exact = ", -8/3, ")")
```

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
\frac{\partial I}{\partial\lambda} = -\frac{2}{(1+\mu)^2(1+\mu+1)},
\quad
\frac{\partial I}{\partial\mu} = \frac{-2}{(1+\mu)^3}
```

```@example shared
f9 = (u, p) -> (u[1] * u[2])^p[2] * exp(-p[1] * (u[1] + u[2]))
prob9 = IntegralProblem(f9, (zeros(2), ones(2)), [0.0, -0.3])

# I as a function of p = [λ, μ] — singularity exponent = p[2], held fixed
Ivec(p) = solve(remake(prob9, p=p),
                DecuhrAlgorithm(singul=2, alpha=-0.3);   # alpha fixed at μ₀
                abstol=1e-7).u

p0     = [0.0, -0.3]
grad   = ForwardDiff.gradient(Ivec, p0)

μ = p0[2]
exact_I    = 1/(1+μ)^2
exact_dIdλ = -2/(1+μ)^2 / (2+μ)     # -(∫ x^{1/2} dx)² = -(2/(2+μ)) · ...
# dI/dλ|λ=0 = -∫∫ (x₁+x₂)(x₁x₂)^μ dx = -2/(μ+2)/(μ+1)
exact_dIdλ = -2 / ((1+μ) * (2+μ))
# dI/dμ = d/dμ [1/(1+μ)²] = -2/(1+μ)³
exact_dIdμ = -2/(1+μ)^3

println("∇I ≈ ", grad)
println("exact ∇I = [", exact_dIdλ, ", ", exact_dIdμ, "]")
println("|err|    = ", abs.(grad .- [exact_dIdλ, exact_dIdμ]))
```

!!! warning "Alpha held fixed during gradient computation"
    When differentiating with respect to ``\mu`` (the singularity exponent),
    `alpha` in `DecuhrAlgorithm` must remain a constant `Float64` —
    it controls the extrapolation rule, not the value of the integral.
    For an exact gradient in ``\mu``, one must evaluate at ``\mu_0``
    and accept that the quadrature error introduces a bias of order
    ``O(\text{abstol})``.

---

## 9 — Budget control (`MaxIters`)

Behaviour when the budget is too tight:

```@example shared
sol9 = solve(prob5,
             DecuhrAlgorithm(singul=3, alpha=-1/3);
             abstol=1e-14,   # very tight tolerance
             maxiters=500)   # very limited budget

println("retcode : ", sol9.retcode)       # MaxIters
println("best estimate : ", sol9.u)      # value still available
println("|err| ≈ ", abs(sol9.u - (3/2)^3))
```

---

## 10 — Diagnostics: `sol.stats` and a conservative `MaxIters`

`sol.stats` reports the number of integrand evaluations and the raw DECUHR
code.  Because DECUHR's error estimator is **deliberately conservative**, a
`MaxIters` return code often accompanies a result that is already accurate to
far better than the requested tolerance — always inspect `sol.u`/`sol.resid`
rather than trusting the return code alone.

```@example shared
sol10 = solve(prob, DecuhrAlgorithm(singul=2, alpha=-0.5); abstol=1e-12)

println("retcode   : ", sol10.retcode)        # may be MaxIters at this tolerance
println("numevals  : ", sol10.stats.numevals) # integrand evaluations
println("ifail     : ", sol10.stats.ifail)    # raw DECUHR code
println("I ≈ ", sol10.u, "   |err| = ", abs(sol10.u - 4.0))  # accurate regardless
```
