# DECALP paths: automatic estimation of the singularity exponent alpha and
# automatic detection of the logarithmic factor (Björstad–Grosse–Dahlquist
# elimination), which only run when the supplied alpha satisfies alpha ≤ -singul.

@testset "2-D vertex singularity, alpha auto-estimated" begin
    f(u, _) = (u[1] * u[2])^(-0.5)
    prob = IntegralProblem(f, (zeros(2), ones(2)))
    sol = solve(
        prob,
        DecuhrAlgorithm(singul = 2);
        abstol = 1.0e-7, reltol = 1.0e-6
    )
    @test isapprox(sol.u, 4.0; rtol = 1.0e-3)
end

@testset "Logarithmic factor auto-detected (BGD elimination)" begin
    # -log(x·y): alpha is left at its default (-2 ≤ -singul), so DECALP must
    # both estimate alpha ≈ 0 and detect logf = 1 on its own.
    f(u, _) = -log(u[1] * u[2])
    prob = IntegralProblem(f, (zeros(2), ones(2)))
    sol = solve(prob, DecuhrAlgorithm(singul = 2); abstol = 1.0e-8, reltol = 1.0e-6)
    @test isapprox(sol.u, 2.0; rtol = 1.0e-5)
end

@testset "Auto-estimated alpha matches the explicit-alpha result" begin
    f(u, _) = (u[1] * u[2])^(-0.5)
    prob = IntegralProblem(f, (zeros(2), ones(2)))
    kwargs = (abstol = 1.0e-7, reltol = 1.0e-6, maxiters = 300_000)
    auto = solve(prob, DecuhrAlgorithm(singul = 2); kwargs...)
    expl = solve(prob, DecuhrAlgorithm(singul = 2, alpha = -0.5); kwargs...)
    # The two paths subdivide differently; they agree to the requested accuracy,
    # not bit-for-bit (measured gap ≈ 1.5e-5 on 4.0).
    @test isapprox(auto.u, expl.u; rtol = 1.0e-4)
end
