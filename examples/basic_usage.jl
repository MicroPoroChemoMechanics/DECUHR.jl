# basic_usage.jl — DECUHR.jl usage examples
#
# Demonstrates the DECUHR algorithm as a pluggable Integrals.jl solver.
# Run from the DECUHR.jl directory with:
#
#   julia --project=. examples/basic_usage.jl
#
# or activate the environment first:
#
#   julia> using Pkg; Pkg.activate("."); Pkg.instantiate()
#   julia> include("examples/basic_usage.jl")

using Pkg; Pkg.activate(joinpath(@__DIR__, "..")); Pkg.instantiate()

using Integrals
using DECUHR

println("=" ^ 60)
println("DECUHR.jl — Integration examples")
println("=" ^ 60)

# ──────────────────────────────────────────────────────────────
# Example 1 : 2-D vertex singularity  (alpha known)
#
#   ∫₀¹∫₀¹ (x·y)^{-1/2} dx dy = (∫₀¹ x^{-1/2} dx)² = 2² = 4
# ──────────────────────────────────────────────────────────────
println("\n--- Example 1: 2-D vertex singularity, alpha = -0.5 ---")

f1(u, _) = (u[1] * u[2])^(-0.5)

prob1 = IntegralProblem(f1, (zeros(2), ones(2)))
sol1  = solve(prob1,
              DecuhrAlgorithm(singul=2, alpha=-0.5);
              abstol=1e-8, reltol=1e-6)

exact1 = 4.0
@show sol1.u
@show sol1.resid
@show abs(sol1.u - exact1)
@show sol1.retcode
println("Error vs exact (4.0): ", abs(sol1.u - exact1))

solve(IntegralProblem((x, _) -> log(x[1])*log(x[2]), (zeros(2), ones(2))), DecuhrAlgorithm())

# ──────────────────────────────────────────────────────────────
# Example 2 : 2-D vertex singularity  (alpha auto-estimated)
#
#   ∫₀¹∫₀¹ (x·y)^{-1/2} dx dy = 4  (same, but alpha auto)
# ──────────────────────────────────────────────────────────────
println("\n--- Example 2: same integral, alpha auto-estimated ---")

sol2 = solve(prob1,
             DecuhrAlgorithm(singul=2);   # alpha=-2 ≤ -singul → auto
             abstol=1e-7, reltol=1e-6)

@show sol2.u
@show abs(sol2.u - exact1)
@show sol2.retcode

# ──────────────────────────────────────────────────────────────
# Example 3 : regular 2-D integrand (no singularity)
#
#   ∫₀^{π/2}∫₀^{π/2} sin(x)·cos(y) dx dy = 1
# ──────────────────────────────────────────────────────────────
println("\n--- Example 3: smooth 2-D integrand ---")

f3(u, _) = sin(u[1]) * cos(u[2])

prob3 = IntegralProblem(f3, (zeros(2), fill(π/2, 2)))
sol3  = solve(prob3, DecuhrAlgorithm(); abstol=1e-10, reltol=1e-10)

exact3 = 1.0
@show sol3.u
@show abs(sol3.u - exact3)
@show sol3.retcode

# ──────────────────────────────────────────────────────────────
# Example 4 : vectorial integration (NUMFUN = 2)
#
#   ∫₀¹∫₀¹ x² + y²  dx dy = 2/3
#   ∫₀¹∫₀¹ x · y    dx dy = 1/4
# ──────────────────────────────────────────────────────────────
println("\n--- Example 4: vector integrand (NUMFUN=2) ---")

f4(u, _) = [u[1]^2 + u[2]^2,  u[1]*u[2]]

prob4 = IntegralProblem(f4, (zeros(2), ones(2)))
sol4  = solve(prob4, DecuhrAlgorithm(); abstol=1e-9)

exact4 = [2/3, 1/4]
@show sol4.u
@show abs.(sol4.u .- exact4)
@show sol4.retcode

# ──────────────────────────────────────────────────────────────
# Example 5 : 3-D weak singularity at origin
#
#   ∫₀¹∫₀¹∫₀¹ (x·y·z)^{-1/3} dx dy dz = (3/2)³ = 27/8
# ──────────────────────────────────────────────────────────────
println("\n--- Example 5: 3-D vertex singularity, alpha = -1/3 ---")

f5(u, _) = (u[1] * u[2] * u[3])^(-1/3)

prob5 = IntegralProblem(f5, (zeros(3), ones(3)))
sol5  = solve(prob5,
              DecuhrAlgorithm(singul=3, alpha=-1/3);
              abstol=1e-7, reltol=1e-7)

exact5 = (3/2)^3   # = 27/8 = 3.375
@show sol5.u
@show abs(sol5.u - exact5)
@show sol5.retcode

# ──────────────────────────────────────────────────────────────
# Example 6 : logarithmic singularity
#
#   ∫₀¹∫₀¹ -log(x·y)  dx dy = 2   (singularity at origin)
#   Note: -log(x·y) = -(log x + log y) = |log x| + |log y| → +∞ at 0
# ──────────────────────────────────────────────────────────────
println("\n--- Example 6: logarithmic singularity ---")

f6(u, _) = -log(u[1] * u[2])   # = -log(x) - log(y)

prob6 = IntegralProblem(f6, (zeros(2), ones(2)))
sol6  = solve(prob6,
              DecuhrAlgorithm(singul=2, alpha=0.0, logf=1);
              abstol=1e-8)

exact6 = 2.0
@show sol6.u
@show abs(sol6.u - exact6)
@show sol6.retcode

# ──────────────────────────────────────────────────────────────
# Example 7 : MaxIters behaviour
# ──────────────────────────────────────────────────────────────
println("\n--- Example 7: maxiters too small → MaxIters retcode ---")

sol7 = solve(prob5,
             DecuhrAlgorithm(singul=3, alpha=-1/3);
             abstol=1e-14,   # very tight
             maxiters=500)   # very tight budget

@show sol7.retcode
@show sol7.u   # still returns best estimate

println("\n" * "=" ^ 60)
println("All examples completed.")
println("=" ^ 60)
