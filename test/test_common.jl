import Ogre
using FactCheck

facts("Common functionality") do
    context("Value sets") do
        vs1 = Ogre.ValueSet(1,2,3)
        vs2 = Ogre.ValueSet(1,2,3,4)
        @fact_throws ValueSet(1,2)

        @fact Ogre.mass(vs1) => 1
        @fact Ogre.nonmass(vs1) => [2,3]
        @fact Ogre.mass(vs2) => 1
        @fact Ogre.nonmass(vs2) => [2,3,4]
    end

    context("Copy-modification of types") do
        vs1 = Ogre.ValueSet(1,2,3)
        vs2 = Ogre.cpmod(vs1, r=4)

        @fact vs2 => Ogre.ValueSet(1,4,3)
    end
end
