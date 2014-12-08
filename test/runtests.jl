using FactCheck

FactCheck.setstyle(:compact)

include("test_eos.jl")
include("test_structure.jl")
include("test_integrator.jl")

using Lint

println("Linting...")
lintpkg("ogre")

FactCheck.exitstatus()
