# Test runner script
using Base.Test

# tests

@testset "Ogre tests" begin
    include("test_common.jl")
    include("test_config.jl")
        include("test_config_defaults.jl")
    include("test_constants.jl")
    include("test_eos.jl")
    include("test_heatcapacity.jl")
    include("test_integrator.jl")
    include("test_physicalvalues.jl")
    include("test_structure.jl")
    include("test_util.jl")
end
