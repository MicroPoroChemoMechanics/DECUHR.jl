# Integrals.jl-facing behaviour: construction, representative integrals of each
# kind (singular, radial, smooth, logarithmic, polynomial), rule keys and stats.

@testset "Smoke: algorithm construction" begin
    @test DecuhrAlgorithm() isa DecuhrAlgorithm
    @test DecuhrAlgorithm(singul = 2, alpha = -0.5) isa DecuhrAlgorithm
end

@testset "2-D vertex singularity, alpha = -0.5" begin
    # ∫₀¹∫₀¹ (x·y)^(-1/2) dx dy = (∫₀¹ x^(-1/2) dx)² = 4
    f(u, _) = (u[1] * u[2])^(-0.5)
    prob = IntegralProblem(f, (zeros(2), ones(2)))
    sol = solve(
        prob,
        DecuhrAlgorithm(singul = 2, alpha = -0.5);
        abstol = 1.0e-8, reltol = 1.0e-6, maxiters = 300_000
    )
    @test isapprox(sol.u, 4.0; rtol = 1.0e-3)
end

@testset "2-D radial vertex singularity (1/√(x²+y²))" begin
    # ∫₀¹∫₀¹ 1/√(x²+y²) dx dy = 2·ln(1+√2)
    f(u, _) = 1 / sqrt(u[1]^2 + u[2]^2)
    prob = IntegralProblem(f, (zeros(2), ones(2)))
    sol = solve(
        prob,
        DecuhrAlgorithm(singul = 2, alpha = -1.0);
        abstol = 1.0e-8, reltol = 1.0e-7
    )
    @test isapprox(sol.u, 2 * log(1 + sqrt(2)); rtol = 1.0e-5)
end

@testset "Smooth 2-D integrand" begin
    # ∫₀^{π/2}∫₀^{π/2} sin(x)·cos(y) dx dy = 1
    f(u, _) = sin(u[1]) * cos(u[2])
    prob = IntegralProblem(f, (zeros(2), fill(π / 2, 2)))
    sol = solve(prob, DecuhrAlgorithm(); abstol = 1.0e-10, reltol = 1.0e-10)
    @test isapprox(sol.u, 1.0; atol = 1.0e-9)
end

@testset "3-D vertex singularity, alpha = -1/3" begin
    # ∫₀¹∫₀¹∫₀¹ (x·y·z)^(-1/3) dx dy dz = (3/2)³ = 27/8
    f(u, _) = (u[1] * u[2] * u[3])^(-1 / 3)
    prob = IntegralProblem(f, (zeros(3), ones(3)))
    sol = solve(
        prob,
        DecuhrAlgorithm(singul = 3, alpha = -1 / 3, wrksub = 60_000);
        abstol = 1.0e-7, reltol = 1.0e-7, maxiters = 1_500_000
    )
    @test isapprox(sol.u, 27 / 8; rtol = 5.0e-3)
end

@testset "Logarithmic singularity" begin
    # ∫₀¹∫₀¹ -log(x·y) dx dy = 2
    f(u, _) = -log(u[1] * u[2])
    prob = IntegralProblem(f, (zeros(2), ones(2)))
    sol = solve(
        prob,
        DecuhrAlgorithm(singul = 2, alpha = 0.0, logf = 1);
        abstol = 1.0e-8
    )
    @test isapprox(sol.u, 2.0; atol = 1.0e-6)
end

@testset "Polynomial integrated exactly by the cubature rule" begin
    # A genuine Genz–Malik rule integrates a degree-2 polynomial exactly.
    # ∫₀¹∫₀¹ (x²+y²) dx dy = 2/3 — should hit machine precision.
    f(u, _) = u[1]^2 + u[2]^2
    prob = IntegralProblem(f, (zeros(2), ones(2)))
    sol = solve(prob, DecuhrAlgorithm(); abstol = 1.0e-12, reltol = 1.0e-12)
    @test isapprox(sol.u, 2 / 3; atol = 1.0e-12)
end

@testset "Rule key sweep on a smooth integrand" begin
    # x²+y² is degree 2 → exact for every rule (keys 1, 3, 4 in 2-D).
    f2(u, _) = u[1]^2 + u[2]^2
    prob2 = IntegralProblem(f2, (zeros(2), ones(2)))
    for key in (1, 3, 4)
        sol = solve(prob2, DecuhrAlgorithm(key = key); abstol = 1.0e-10)
        @test isapprox(sol.u, 2 / 3; atol = 1.0e-9)
    end
    # key = 2 is the 3-D (degree-11) rule.
    f3(u, _) = u[1]^2 + u[2]^2 + u[3]^2
    prob3 = IntegralProblem(f3, (zeros(3), ones(3)))
    sol3 = solve(prob3, DecuhrAlgorithm(key = 2); abstol = 1.0e-10)
    @test isapprox(sol3.u, 1.0; atol = 1.0e-9)   # ∫(x²+y²+z²) over unit cube = 1
end

@testset "Solution stats expose neval and ifail" begin
    f(u, _) = sin(u[1]) * cos(u[2])
    prob = IntegralProblem(f, (zeros(2), fill(π / 2, 2)))
    sol = solve(prob, DecuhrAlgorithm(); abstol = 1.0e-9)
    @test sol.stats.numevals > 0
    @test sol.stats.ifail == 0
    # Singular case that exhausts the default budget keeps a usable result.
    fs(u, _) = (u[1] * u[2])^(-0.5)
    probs = IntegralProblem(fs, (zeros(2), ones(2)))
    sols = solve(probs, DecuhrAlgorithm(singul = 2, alpha = -0.5); abstol = 1.0e-12)
    @test sols.stats.numevals > 0
    @test isapprox(sols.u, 4.0; rtol = 1.0e-3)   # correct even if retcode == MaxIters
end

@testset "Extrapolation weighted-correction branch (small emax)" begin
    # emax = 3 forces NUMU > EMAX on a strongly singular integrand, which
    # exercises the weighted-correction branch (and its BETA-index clamp).
    f(u, _) = (u[1] * u[2])^(-0.5)
    prob = IntegralProblem(f, (zeros(2), ones(2)))
    sol = solve(
        prob,
        DecuhrAlgorithm(singul = 2, alpha = -0.5, emax = 3);
        abstol = 1.0e-6, reltol = 1.0e-5, maxiters = 200_000
    )
    @test isapprox(sol.u, 4.0; rtol = 2.0e-2)
end
