using FactCheck
FactCheck.setstyle(:default) # :compact or :default
FactCheck.clear_results()

function main()
    include("test_eos.jl")
    include("test_structure.jl")
    include("test_integrator.jl")

    FactCheck.exitstatus()
end

main()

