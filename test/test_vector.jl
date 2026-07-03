# Vector integrands (NUMFUN ≥ 2): all components share the same subdivision.

@testset "Vector integrand (NUMFUN = 2)" begin
    # ∫₀¹∫₀¹ [x²+y², x·y] = [2/3, 1/4]
    f(u, _) = [u[1]^2 + u[2]^2, u[1] * u[2]]
    prob = IntegralProblem(f, (zeros(2), ones(2)))
    sol = solve(prob, DecuhrAlgorithm(); abstol = 1.0e-9)
    @test isapprox(sol.u[1], 2 / 3; atol = 1.0e-7)
    @test isapprox(sol.u[2], 1 / 4; atol = 1.0e-7)
end

@testset "In-place vector integrand (IntegralFunction{true})" begin
    # Same integrals through the in-place SciML form: f(y, u, p) fills y.
    # numfun and the value type come from the prototype — no probe evaluation.
    f!(y, u, _) = (y[1] = u[1]^2 + u[2]^2; y[2] = u[1] * u[2]; nothing)
    fun = IntegralFunction(f!, zeros(2))
    prob = IntegralProblem(fun, (zeros(2), ones(2)))
    sol = solve(prob, DecuhrAlgorithm(); abstol = 1.0e-9)
    @test isapprox(sol.u[1], 2 / 3; atol = 1.0e-7)
    @test isapprox(sol.u[2], 1 / 4; atol = 1.0e-7)

    # The in-place path is the zero-allocation route: a second solve must stay
    # well below the per-evaluation cost of the out-of-place form (one fresh
    # vector per call).  Bound total allocated bytes by a generous constant
    # rather than zero to stay robust across Julia versions and CI settings.
    solve(prob, DecuhrAlgorithm(); abstol = 1.0e-9)   # warm-up
    fun_oop = IntegralFunction((u, _) -> [u[1]^2 + u[2]^2, u[1] * u[2]])
    prob_oop = IntegralProblem(fun_oop, (zeros(2), ones(2)))
    solve(prob_oop, DecuhrAlgorithm(); abstol = 1.0e-9)   # warm-up
    bytes_iip = @allocated solve(prob, DecuhrAlgorithm(); abstol = 1.0e-9)
    bytes_oop = @allocated solve(prob_oop, DecuhrAlgorithm(); abstol = 1.0e-9)
    @test bytes_iip < bytes_oop
end

@testset "Singular in-place vector integrand" begin
    f!(y, u, _) = (y[1] = (u[1] * u[2])^(-0.5); y[2] = u[1] + u[2]; nothing)
    fun = IntegralFunction(f!, zeros(2))
    prob = IntegralProblem(fun, (zeros(2), ones(2)))
    sol = solve(
        prob,
        DecuhrAlgorithm(singul = 2, alpha = -0.5);
        abstol = 1.0e-8, reltol = 1.0e-6, maxiters = 1_000_000
    )
    @test isapprox(sol.u[1], 4.0; rtol = 1.0e-3)
    @test isapprox(sol.u[2], 1.0; rtol = 1.0e-6)
end

@testset "Singular vector integrand (NUMFUN = 2)" begin
    # A singular and a smooth component side by side:
    # ∫₀¹∫₀¹ [(x·y)^(-1/2), x+y] = [4, 1].
    f(u, _) = [(u[1] * u[2])^(-0.5), u[1] + u[2]]
    prob = IntegralProblem(f, (zeros(2), ones(2)))
    sol = solve(
        prob,
        DecuhrAlgorithm(singul = 2, alpha = -0.5);
        abstol = 1.0e-8, reltol = 1.0e-6, maxiters = 1_000_000
    )
    @test sol.u isa AbstractVector
    @test length(sol.u) == 2
    @test isapprox(sol.u[1], 4.0; rtol = 1.0e-3)
    @test isapprox(sol.u[2], 1.0; rtol = 1.0e-6)
end
