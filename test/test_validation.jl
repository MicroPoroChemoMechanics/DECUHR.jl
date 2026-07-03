# Parameter validation: every reachable DECHHR failure code. The raw code is
# exposed as sol.stats.ifail and mapped to ReturnCode.Failure.

@testset "Parameter validation → ReturnCode.Failure" begin
    f(u, _) = 1.0
    prob2 = IntegralProblem(f, (zeros(2), ones(2)))
    prob3 = IntegralProblem(f, (zeros(3), ones(3)))

    # ifail = 2 — KEY out of range
    sol = solve(prob2, DecuhrAlgorithm(key = 5))
    @test sol.retcode == RC.Failure
    @test sol.stats.ifail == 2
    @test sol.stats.message == DECUHR.ifail_message(2)
    # ifail = 3 — NDIM out of range (1 and > 15)
    @test solve(IntegralProblem(f, ([0.0], [1.0])), DecuhrAlgorithm()).stats.ifail == 3
    @test solve(IntegralProblem(f, (zeros(16), ones(16))), DecuhrAlgorithm()).stats.ifail == 3
    # ifail = 4 — KEY = 1 requires NDIM = 2
    @test solve(prob3, DecuhrAlgorithm(key = 1)).stats.ifail == 4
    # ifail = 5 — KEY = 2 requires NDIM = 3
    @test solve(prob2, DecuhrAlgorithm(key = 2)).stats.ifail == 5
    # ifail = 6 — NUMFUN < 1 (empty integrand output)
    @test solve(IntegralProblem((u, _) -> Float64[], (zeros(2), ones(2))), DecuhrAlgorithm()).stats.ifail == 6
    # ifail = 8 — MAXPTS < MINPTS
    @test solve(prob2, DecuhrAlgorithm(minpts = 1000); maxiters = 100).stats.ifail == 8
    # ifail = 9 — both tolerances negative
    @test solve(prob2, DecuhrAlgorithm(); abstol = -1.0, reltol = -1.0).stats.ifail == 9
    # ifail = 11 — SINGUL > NDIM
    @test solve(prob2, DecuhrAlgorithm(singul = 3)).stats.ifail == 11
    # ifail = 12 — LOGF outside {0, 1}
    @test solve(prob2, DecuhrAlgorithm(singul = 2, alpha = -0.5, logf = 2)).stats.ifail == 12
    # ifail = 13 — MAXPTS too small for even the first rule applications
    @test solve(prob2, DecuhrAlgorithm(); maxiters = 10).stats.ifail == 13
    # ifail = 15 — EMAX < 1
    @test solve(prob2, DecuhrAlgorithm(emax = 0)).stats.ifail == 15
    # ifail = 17 — WRKSUB too small for the requested budget
    @test solve(prob2, DecuhrAlgorithm(wrksub = 10)).stats.ifail == 17
end

@testset "Scalar domain raises a clear error" begin
    # DECUHR is 2 ≤ ndim ≤ 15 only; a scalar domain used to die with an obscure
    # MethodError deep in the driver.
    f(u, _) = 1.0
    prob = IntegralProblem(f, (0.0, 1.0))
    @test_throws ArgumentError solve(prob, DecuhrAlgorithm())
end

@testset "ifail = 7 (bad bounds) at the driver level" begin
    # Integrals.jl's change of variables absorbs reversed bounds before DECUHR
    # sees them, so this code is only reachable through the driver.
    funsub = (x, fv) -> (fv[1] = 1.0; nothing)
    _, _, _, ifail = DECUHR._decuhr_driver(
        2, 1, [1.0, 1.0], [0.0, 0.0], 0, 100_000, funsub, Float64,
        1, 0.0, 0, 1.0e-8, 1.0e-6, 0, 50_000, 20
    )
    @test ifail == 7
end

@testset "ifail_message covers every reachable code" begin
    for code in (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 17)
        @test DECUHR.ifail_message(code) isa String
        @test !startswith(DECUHR.ifail_message(code), "unknown")
    end
    @test startswith(DECUHR.ifail_message(99), "unknown")
end

@testset "Infinite domains keep the standard transformation path" begin
    # The finite-domain ChangeOfVariables bypass must not intercept infinite
    # bounds; those keep Integrals' remap and fail exactly as they always did
    # (DECUHR's vertex-singularity model does not apply to infinite domains).
    f(u, _) = exp(-u[1]^2 - u[2]^2)
    prob = IntegralProblem(f, ([-Inf, -Inf], [Inf, Inf]))
    sol = solve(prob, DecuhrAlgorithm(); abstol = 1.0e-8, reltol = 1.0e-6)
    @test sol.retcode == RC.Failure
    @test sol.stats.ifail == 14
end

@testset "ifail = 14 — alpha estimation failure" begin
    # An identically zero integrand defeats DECALP (no usable ratio table).
    sol = solve(IntegralProblem((u, _) -> 0.0, (zeros(2), ones(2))), DecuhrAlgorithm(singul = 2))
    @test sol.retcode == RC.Failure
    @test sol.stats.ifail == 14
end
