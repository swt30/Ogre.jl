using FactCheck
using Lint

facts("OGRE package tests") do
    include("test_eos.jl")
    include("test_structure.jl")
    include("test_integrator.jl")
end
