using FactCheck
FactCheck.setstyle(:compact) # :compact or :default

function main()
    include("test_eos.jl")
    include("test_structure.jl")
    include("test_integrator.jl")

    FactCheck.exitstatus()
end

main()

