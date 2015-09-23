# Test runner script

using FactCheck
    FactCheck.setstyle(:default)
    FactCheck.clear_results()


# tests

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

FactCheck.exitstatus()
