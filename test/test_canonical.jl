# Canonical validation cases C1–C6 (Espelid & Genz 1994), locked BIT-FOR-BIT.
#
# The golden values below were captured from DECUHR.jl v0.1.3 on x86-64 /
# Julia 1.12 (openlibm) with the exact parameters of fortran_ref/COMPARAISON.md
# (abstol = 1e-8, reltol = 1e-6, maxpts = 1e6, wrksub = 50000, emax = 20,
# automatic rule). At the driver level they agree with the original Fortran
# DECUHR (5/6 cases to the last bit; C4 differs by ~1.3e-10 from Fortran due to
# a documented floating-point reassociation — the value locked here is DECUHR.jl's own).
#
# Results are compared by reinterpreting the Float64 bits: any change in the
# numerical path (reordered accumulation, @simd/@fastmath, different libm,
# altered subdivision) must fail here and be a deliberate, documented decision.
#
# Two levels are locked:
#   * driver — DECUHR._decuhr_driver, the pure algorithm (the Fortran contract);
#   * solve  — the user-facing Integrals.jl path.  Since DECUHR skips the
#     Integrals ChangeOfVariables remap on finite domains, solve agrees with
#     the driver bit-for-bit (same result, abserr, neval, ifail).

@testset "Canonical C1-C6 (bit-exact golden)" begin
    maxpts = 1_000_000
    wrksub = 50_000
    emax = 20
    abstol = 1.0e-8
    reltol = 1.0e-6

    cases = [
        (
            name = "C1 (x·y)^(-1/2)",
            f = u -> (u[1] * u[2])^(-0.5),
            lb = [0.0, 0.0], ub = [1.0, 1.0], singul = 2, alpha = -0.5, logf = 0,
            driver = (
                result = 0x4010000feece1bc9,   # 4.0000607789324443
                abserr = 0x3ed0835aa7c6d0a9,   # 3.9370303252414573e-6
                neval = 110045, ifail = 0,
            ),
            solve = (result = 0x4010000feece1bc9, neval = 110045, ifail = 0),
        ),
        (
            name = "C2 1/√(x²+y²)",
            f = u -> 1 / sqrt(u[1]^2 + u[2]^2),
            lb = [0.0, 0.0], ub = [1.0, 1.0], singul = 2, alpha = -1.0, logf = 0,
            driver = (
                result = 0x3ffc34365f568eb9,   # 1.7627471660752418
                abserr = 0x3e61a7ef1c950400,   # 3.2886848420298463e-8
                neval = 715, ifail = 0,
            ),
            solve = (result = 0x3ffc34365f568eb9, neval = 715, ifail = 0),
        ),
        (
            name = "C3 sin(x)·cos(y)",
            f = u -> sin(u[1]) * cos(u[2]),
            lb = [0.0, 0.0], ub = [pi / 2, pi / 2], singul = 1, alpha = 0.0, logf = 0,
            driver = (
                result = 0x3ff0000000000034,   # 1.0000000000000115
                abserr = 0x3d4fb9d6ec58682c,   # 2.2542633159715836e-13
                neval = 195, ifail = 0,
            ),
            solve = (result = 0x3ff0000000000034, neval = 195, ifail = 0),
        ),
        (
            name = "C4 (x·y·z)^(-1/3)",
            f = u -> (u[1] * u[2] * u[3])^(-1 / 3),
            lb = [0.0, 0.0, 0.0], ub = [1.0, 1.0, 1.0], singul = 3, alpha = -1 / 3, logf = 0,
            driver = (
                result = 0x400b00103fef890b,   # 3.3750309939360981
                abserr = 0x3f3a31ce04184351,   # 3.9969711005576498e-4
                neval = 999871, ifail = 1,
            ),
            solve = (result = 0x400b00103fef890b, neval = 999871, ifail = 1),
        ),
        (
            name = "C5 x²+y²",
            f = u -> u[1]^2 + u[2]^2,
            lb = [0.0, 0.0], ub = [1.0, 1.0], singul = 1, alpha = 0.0, logf = 0,
            driver = (
                result = 0x3fe5555555555556,   # 0.66666666666666674
                abserr = 0x3d005496840f2373,   # 7.2521640646959477e-15
                neval = 195, ifail = 0,
            ),
            solve = (result = 0x3fe5555555555556, neval = 195, ifail = 0),
        ),
        (
            name = "C6 -log(x·y)",
            f = u -> -log(u[1] * u[2]),
            lb = [0.0, 0.0], ub = [1.0, 1.0], singul = 2, alpha = 0.0, logf = 1,
            driver = (
                result = 0x4000000003dd3e52,   # 2.0000000287907733
                abserr = 0x3ebfde3c44736c2e,   # 1.899487203796715e-6
                neval = 9035, ifail = 0,
            ),
            solve = (result = 0x4000000003dd3e52, neval = 9035, ifail = 0),
        ),
    ]

    @testset "driver: $(c.name)" for c in cases
        funsub = (x, fv) -> (fv[1] = c.f(x); nothing)
        r, e, neval, ifail = DECUHR._decuhr_driver(
            length(c.lb), 1, c.lb, c.ub, 0, maxpts, funsub, Float64,
            c.singul, c.alpha, c.logf, abstol, reltol, 0, wrksub, emax
        )
        @test reinterpret(UInt64, r[1]) == c.driver.result
        @test reinterpret(UInt64, e[1]) == c.driver.abserr
        @test neval == c.driver.neval
        @test ifail == c.driver.ifail
    end

    @testset "solve: $(c.name)" for c in cases
        f = (u, _) -> c.f(u)
        prob = IntegralProblem(f, (c.lb, c.ub))
        alg = DecuhrAlgorithm(
            singul = c.singul, alpha = c.alpha, logf = c.logf,
            wrksub = wrksub, emax = emax
        )
        sol = solve(prob, alg; abstol = abstol, reltol = reltol, maxiters = maxpts)
        @test reinterpret(UInt64, sol.u) == c.solve.result
        @test sol.stats.numevals == c.solve.neval
        @test sol.stats.ifail == c.solve.ifail
        @test sol.retcode == (c.solve.ifail == 0 ? RC.Success : RC.MaxIters)
    end
end
