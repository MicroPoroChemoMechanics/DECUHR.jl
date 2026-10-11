# Automatic differentiation with respect to the problem parameters p, with and
# without automatic alpha estimation (which always runs on primal values).

@testset "ForwardDiff: derivative w.r.t. a parameter (explicit alpha)" begin
    # ∫₀¹∫₀¹ p·(x·y)^(-1/2) dx dy = 4p ⇒ d/dp = 4
    g = ForwardDiff.derivative(2.0) do p
        f = (u, _) -> p * (u[1] * u[2])^(-0.5)
        prob = IntegralProblem(f, (zeros(2), ones(2)))
        sol = solve(
            prob, DecuhrAlgorithm(singul = 2, alpha = -0.5);
            abstol = 1.0e-8, reltol = 1.0e-6, maxiters = 300_000
        )
        sol.u
    end
    @test isapprox(g, 4.0; rtol = 1.0e-2)
end

@testset "ForwardDiff: derivative with auto-estimated alpha" begin
    # Same integral, but alpha is NOT supplied: it must be auto-estimated on
    # the primal integrand, and solve must remain differentiable.
    # ∫₀¹∫₀¹ p·(x·y)^(-1/2) = 4p ⇒ d/dp = 4.
    g = ForwardDiff.derivative(2.0) do p
        f = (u, _) -> p * (u[1] * u[2])^(-0.5)
        prob = IntegralProblem(f, (zeros(2), ones(2)))
        sol = solve(
            prob, DecuhrAlgorithm(singul = 2);
            abstol = 1.0e-7, reltol = 1.0e-6, maxiters = 300_000
        )
        sol.u
    end
    @test isapprox(g, 4.0; rtol = 1.0e-2)
end

@testset "ForwardDiff: nested derivatives with respect to a parameter" begin
    # ∫₀¹∫₀¹ (x·y)^(-1/2) e^(p·x) dx dy = 2 ∫₀¹ x^(-1/2) e^(p·x) dx, so at p = 0
    # its first and second derivatives are 2∫x^(1/2) = 4/3 and 2∫x^(3/2) = 4/5.
    I(p) = solve(
        IntegralProblem((u, q) -> (u[1] * u[2])^(-0.5) * exp(q * u[1]), (zeros(2), ones(2)), p),
        DecuhrAlgorithm(singul = 2, alpha = -0.5);
        abstol = 1.0e-10, reltol = 1.0e-10, maxiters = 2_000_000
    ).u
    @test ForwardDiff.derivative(I, 0.0) ≈ 4 / 3 rtol = 1.0e-10
    @test ForwardDiff.derivative(p -> ForwardDiff.derivative(I, p), 0.0) ≈ 4 / 5 rtol = 1.0e-10
end

@testset "ForwardDiff: the bounds are refused, their mapping is the way" begin
    # ∫₀ᵇ∫₀¹ (x·y)^(-1/2) = 4√b, d/db = 2/√b. A dual bound is refused by name;
    # mapped onto the unit square, x = b·t with Jacobian b, b is a parameter.
    f(u, _) = (u[1] * u[2])^(-0.5)
    alg = DecuhrAlgorithm(singul = 2, alpha = -0.5)
    kw = (abstol = 1.0e-10, reltol = 1.0e-10, maxiters = 2_000_000)
    J(b) = solve(IntegralProblem(f, ([0.0, 0.0], [b, 1.0])), alg; kw...).u
    @test_throws "does not differentiate with respect to the bounds" ForwardDiff.derivative(J, 1.0)
    # The derivative of the mapped integrand, √b·(t·y)^(-1/2), is the integrand
    # itself over 2b: it carries the accuracy of the value, 6e-8 relative in
    # this budget, not that of the nested-derivative test above.
    Jmapped(b) = solve(IntegralProblem((t, b) -> b * f([b * t[1], t[2]], nothing), (zeros(2), ones(2)), b), alg; kw...).u
    @test ForwardDiff.derivative(Jmapped, 1.0) ≈ 2 rtol = 1.0e-6
end
