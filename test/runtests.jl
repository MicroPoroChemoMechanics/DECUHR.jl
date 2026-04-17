using DECUHR
using Integrals
using Test

@testset "DECUHR.jl" begin

    @testset "Smoke: algorithm construction" begin
        @test DecuhrAlgorithm() isa DecuhrAlgorithm
        @test DecuhrAlgorithm(singul = 2, alpha = -0.5) isa DecuhrAlgorithm
    end

    @testset "2-D vertex singularity, alpha = -0.5" begin
        # ∫₀¹∫₀¹ (x·y)^(-1/2) dx dy = (∫₀¹ x^(-1/2) dx)² = 4
        f(u, _) = (u[1] * u[2])^(-0.5)
        prob = IntegralProblem(f, zeros(2), ones(2))
        sol = solve(
            prob,
            DecuhrAlgorithm(singul = 2, alpha = -0.5);
            abstol = 1.0e-8, reltol = 1.0e-6
        )
        @test isapprox(sol.u, 4.0; rtol = 1.0e-3)
    end

    @testset "2-D vertex singularity, alpha auto-estimated" begin
        f(u, _) = (u[1] * u[2])^(-0.5)
        prob = IntegralProblem(f, zeros(2), ones(2))
        sol = solve(
            prob,
            DecuhrAlgorithm(singul = 2);
            abstol = 1.0e-7, reltol = 1.0e-6
        )
        @test isapprox(sol.u, 4.0; rtol = 1.0e-3)
    end

    @testset "Smooth 2-D integrand" begin
        # ∫₀^{π/2}∫₀^{π/2} sin(x)·cos(y) dx dy = 1
        f(u, _) = sin(u[1]) * cos(u[2])
        prob = IntegralProblem(f, zeros(2), fill(π / 2, 2))
        sol = solve(prob, DecuhrAlgorithm(); abstol = 1.0e-10, reltol = 1.0e-10)
        @test isapprox(sol.u, 1.0; atol = 1.0e-9)
    end

    @testset "3-D vertex singularity, alpha = -1/3" begin
        # ∫₀¹∫₀¹∫₀¹ (x·y·z)^(-1/3) dx dy dz = (3/2)³ = 27/8
        f(u, _) = (u[1] * u[2] * u[3])^(-1 / 3)
        prob = IntegralProblem(f, zeros(3), ones(3))
        sol = solve(
            prob,
            DecuhrAlgorithm(singul = 3, alpha = -1 / 3);
            abstol = 1.0e-7, reltol = 1.0e-7
        )
        @test isapprox(sol.u, 27 / 8; rtol = 5.0e-3)
    end

    @testset "Vector integrand (NUMFUN = 2)" begin
        # ∫₀¹∫₀¹ [x²+y², x·y] = [2/3, 1/4]
        f(u, _) = [u[1]^2 + u[2]^2, u[1] * u[2]]
        prob = IntegralProblem(f, zeros(2), ones(2))
        sol = solve(prob, DecuhrAlgorithm(); abstol = 1.0e-9)
        @test isapprox(sol.u[1], 2 / 3; atol = 1.0e-7)
        @test isapprox(sol.u[2], 1 / 4; atol = 1.0e-7)
    end

    @testset "Logarithmic singularity" begin
        # ∫₀¹∫₀¹ -log(x·y) dx dy = 2
        f(u, _) = -log(u[1] * u[2])
        prob = IntegralProblem(f, zeros(2), ones(2))
        sol = solve(
            prob,
            DecuhrAlgorithm(singul = 2, alpha = 0.0, logf = 1);
            abstol = 1.0e-8
        )
        @test isapprox(sol.u, 2.0; atol = 1.0e-6)
    end

end
