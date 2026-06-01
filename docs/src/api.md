# API Reference

## Module documentation

```@docs
DECUHR.DECUHR
```

## Algorithm type

```@docs
DecuhrAlgorithm
```

## Integrals.jl interface

`DecuhrAlgorithm` integrates with the standard SciML `solve` / `init` / `solve!`
workflow.  No additional API is needed beyond what Integrals.jl provides.

```julia
using Integrals, DECUHR

prob = IntegralProblem(f, (lb, ub))
sol  = solve(prob, DecuhrAlgorithm(; kwargs...);
             abstol=1e-8, reltol=1e-6, maxiters=100_000)
```

### Keyword arguments to `solve`

| Keyword | Default | Description |
| :-------- | :-------- | :------------ |
| `abstol` | `1e-8` | Absolute error tolerance |
| `reltol` | `1e-6` | Relative error tolerance |
| `maxiters` | `100_000` | Maximum number of integrand evaluations |

### Solution fields

| Field | Type | Description |
| :------ | :----- | :------------ |
| `sol.u` | scalar or `Vector` | Integral estimate |
| `sol.resid` | scalar or `Vector` | Absolute error estimate |
| `sol.retcode` | `ReturnCode` | `Success`, `MaxIters`, or `Failure` |

## Usage examples

See the dedicated [Examples](@ref) page for complete, executable examples covering:
singular integrands, vector-valued integrands, parameterised problems,
automatic differentiation with ForwardDiff, and budget-limit behaviour.

## Internal driver (advanced)

```@docs
DECUHR._decuhr_driver
```

!!! note
    `_decuhr_driver` is an internal function exposed for testing and debugging.
    Prefer the standard `solve` interface for everyday use.
