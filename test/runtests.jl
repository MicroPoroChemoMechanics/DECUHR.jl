using Decuhr
using Integrals
using Test
import SciMLBase
import ForwardDiff

const RC = SciMLBase.ReturnCode

@testset "Decuhr.jl" begin
    include("test_canonical.jl")
    include("test_interface.jl")
    include("test_validation.jl")
    include("test_alpha_auto.jl")
    include("test_highdim.jl")
    include("test_vector.jl")
    include("test_forwarddiff.jl")
end
