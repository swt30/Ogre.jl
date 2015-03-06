include("header.jl")

facts("Common functionality") do
    context("Value sets") do
        vs1 = Ogre.ValueSet(1,2,3)
        vs2 = Ogre.ValueSet(1,2,3,4)
        @fact_throws ValueSet(1,2)

        @fact Ogre.mass(vs1) => 1
        @fact Ogre.nonmass(vs1) => [2,3]
        @fact Ogre.mass(vs2) => 1
        @fact Ogre.nonmass(vs2) => [2,3,4]
        @fact Ogre.radius(vs1) => 2
        @fact Ogre.pressure(vs2) => 3
        @fact_throws Ogre.temperature(vs1)
        @fact Ogre.temperature(vs2) => 4
    end

    context("Copy-modification of types") do
        vs1 = Ogre.ValueSet(1,2,3)
        vs1m = Ogre.cpmod(vs1, P=4)
        vs2 = Ogre.ValueSet(1,2,3,4)
        vs2m = Ogre.cpmod(vs2, T=5)

        @fact vs1m => Ogre.ValueSet(1,2,4)
        @fact vs2m => Ogre.ValueSet(1,2,3,5)
    end
end
