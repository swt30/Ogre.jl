include("header.jl")

facts("Common functionality") do
    context("Value sets") do
        vs1 = res.value_set_no_temp
        vs2 = res.value_set_full
        @fact_throws ValueSet(1,2)

        context("Type hierarchy") do
            tempdep = Ogre.ValueSet{Ogre.WithTemp}
            notemp = Ogre.ValueSet{Ogre.NoTemp}
            @fact isa(vs1, tempdep) => false
            @fact isa(vs1, notemp) => true
            @fact isa(vs1, Ogre.ValueSet) => true
            @fact isa(vs2, tempdep) => true
            @fact isa(vs2, notemp) => false
            @fact isa(vs2, Ogre.ValueSet) => true
        end

        context("Get values from them") do
            @fact Ogre.mass(vs1) => 1
            @fact Ogre.nonmass(vs1) => [2,3]
            @fact Ogre.mass(vs2) => 1
            @fact Ogre.nonmass(vs2) => [2,3,4]
            @fact Ogre.radius(vs1) => 2
            @fact Ogre.pressure(vs2) => 3
            @fact_throws Ogre.temperature(vs1)
            @fact Ogre.temperature(vs2) => 4
        end

        context("Zero values or less are considered 'unphysical'") do
            @fact Ogre.zero(Ogre.ValueSet{Ogre.NoTemp}) => not(Ogre.isphysical)
            @fact Ogre.zero(Ogre.ValueSet{Ogre.WithTemp}) => not(Ogre.isphysical)
            @fact Ogre.ValueSet(0, 0, -1) => not(Ogre.isphysical)
            @fact Ogre.ValueSet(0, 0, 0, -1) => not(Ogre.isphysical)
        end

        
    end

    context("Copy-modification of types") do
        vs1 = res.value_set_no_temp
        vs1m = Ogre.cpmod(vs1, P=4)
        vs2 = res.value_set_full
        vs2m = Ogre.cpmod(vs2, T=5)

        @fact vs1m => Ogre.ValueSet(1,2,4)
        @fact vs2m => Ogre.ValueSet(1,2,3,5)
        @pending Ogre.cpmod(vs2, a_field_that_doesnt_exist=6) => "throws"
    end
end
