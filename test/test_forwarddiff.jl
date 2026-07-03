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
