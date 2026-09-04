using DECUHR
using Aqua
using Integrals
using Test
import SciMLBase
import ForwardDiff

const RC = SciMLBase.ReturnCode

"""
    aqua_persistent_tasks(m::Module) -> Bool

`Aqua.test_persistent_tasks` in a child process with default bounds checking.

The check precompiles a synthetic package that depends on `m` and waits for a
sentinel written from that package's module body. Julia writes no
precompilation cache when `--check-bounds` is forced, so the body never runs,
the sentinel never appears, and the check fails for a reason that says nothing
about `m`. `Pkg.test()` forces the flag by default.

Running it in a child with `--check-bounds=auto` keeps both halves: the test
suite proper still runs with bounds checking on, and the check still runs.
"""
function aqua_persistent_tasks(m::Module)
    code = """
    using Aqua, $(nameof(m)), Test
    @testset "persistent_tasks" begin
        Aqua.test_persistent_tasks($(nameof(m)))
    end
    """
    cmd = `$(first(Base.julia_cmd())) --check-bounds=auto --startup-file=no --project=$(Base.active_project()) -e $code`
    return success(run(ignorestatus(cmd)))
end

@testset "DECUHR.jl" begin
    # Ambiguities, unbound type parameters, undefined exports, dependency
    # hygiene and type piracy — none of which the tests below would notice.
    # Nothing is exempted: the package passes every default check as it stands,
    # including the piracy one, whose only candidate here is
    # `Integrals.__solve(..., ::ChangeOfVariables{T, DecuhrAlgorithm}, ...)` —
    # accepted because Aqua descends into type parameters when deciding
    # ownership.
    @testset "Aqua" begin
        # `persistent_tasks` is run apart, see `aqua_persistent_tasks` above.
        Aqua.test_all(DECUHR; persistent_tasks = false)
        @test aqua_persistent_tasks(DECUHR)
    end
    include("test_canonical.jl")
    include("test_interface.jl")
    include("test_validation.jl")
    include("test_alpha_auto.jl")
    include("test_highdim.jl")
    include("test_vector.jl")
    include("test_forwarddiff.jl")
end
