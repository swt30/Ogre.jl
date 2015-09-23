using FactCheck
import Ogre


facts("Value sets") do
    vs_notemp = Ogre.ValueSet(1,2,3)
    vs_full = Ogre.ValueSet(1,2,3,4)
    @fact_throws Ogre.ValueSet(1,2)

    context("Type hierarchy") do
        tempdep = Ogre.ValueSet{Ogre.WithTemp}
        notemp = Ogre.ValueSet{Ogre.NoTemp}
        @fact isa(vs_notemp, tempdep) --> false
        @fact isa(vs_notemp, notemp) --> true
        @fact isa(vs_notemp, Ogre.ValueSet) --> true
        @fact isa(vs_full, tempdep) --> true
        @fact isa(vs_full, notemp) --> false
        @fact isa(vs_full, Ogre.ValueSet) --> true
    end

    context("Get values from them") do
        @fact Ogre.mass(vs_notemp) --> 1
        @fact Ogre.nonmass(vs_notemp) --> [2,3]
        @fact Ogre.mass(vs_full) --> 1
        @fact Ogre.nonmass(vs_full) --> [2,3,4]
        @fact Ogre.radius(vs_notemp) --> 2
        @fact Ogre.pressure(vs_full) --> 3
        @fact_throws Ogre.temperature(vs_notemp)
        @fact Ogre.temperature(vs_full) --> 4
    end

    context("Zero values or less are considered 'unphysical'") do
        @fact Ogre.zero(Ogre.ValueSet{Ogre.NoTemp}) --> not(Ogre.isphysical)
        @fact Ogre.zero(Ogre.ValueSet{Ogre.WithTemp}) --> not(Ogre.isphysical)
        @fact Ogre.ValueSet(0, 0, -1) --> not(Ogre.isphysical)
        @fact Ogre.ValueSet(0, 0, 0, -1) --> not(Ogre.isphysical)
    end
end
