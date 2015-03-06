include("header.jl")

using FactCheck
FactCheck.setstyle(:default) # :compact or :default
FactCheck.clear_results()

function main()

    include("test_common.jl")
    include("test_eos.jl")
    include("test_structure.jl")
    include("test_integrator.jl")
    include("test_heatcapacity.jl")

    FactCheck.exitstatus()
end

main()
