# Dimensions ≥ 4 exercise the general n-D Genz–Malik rules (keys 3 and 4),
# which the 2-D/3-D cases never reach.

@testset "4-D polynomial exactness (keys 3 and 4)" begin
    # ∑xᵢ² has degree 2 → exact for the degree-9 and degree-7 rules.
    f(u, _) = u[1]^2 + u[2]^2 + u[3]^2 + u[4]^2
    prob = IntegralProblem(f, (zeros(4), ones(4)))
    for key in (0, 3, 4)   # key = 0 selects key 3 for ndim ≥ 4
        sol = solve(prob, DecuhrAlgorithm(key = key); abstol = 1.0e-10)
        @test isapprox(sol.u, 4 / 3; atol = 1.0e-9)
    end
end

@testset "4-D edge singularity" begin
    # ∫ x₁^(-1/2) over [0,1]⁴ = 2 (singular on a 3-D face through the vertex).
    f(u, _) = u[1]^(-0.5)
    prob = IntegralProblem(f, (zeros(4), ones(4)))
    sol = solve(
        prob,
        DecuhrAlgorithm(singul = 1, alpha = -0.5);
        abstol = 1.0e-7, reltol = 1.0e-6, maxiters = 1_000_000
    )
    @test isapprox(sol.u, 2.0; rtol = 1.0e-4)
end

@testset "4-D face singularity" begin
    # ∫ (x₁·x₂)^(-1/2) over [0,1]⁴ = 4.
    f(u, _) = (u[1] * u[2])^(-0.5)
    prob = IntegralProblem(f, (zeros(4), ones(4)))
    sol = solve(
        prob,
        DecuhrAlgorithm(singul = 2, alpha = -1.0);
        abstol = 1.0e-7, reltol = 1.0e-6, maxiters = 1_000_000
    )
    @test isapprox(sol.u, 4.0; rtol = 1.0e-3)
end
