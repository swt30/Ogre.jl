include("header.jl")

facts("Heat capacity handling") do
    capfunc(T) = exp(T)

    Cₚ = Ogre.HeatCapacity(Ogre.withtemp, capfunc)
    Cₚ₂ = Ogre.HeatCapacity(capfunc)

    context("Construction") do
        @fact_throws Ogre.HeatCapacity(Ogre.notemp, capfunc)
    end

    context("Values are as expected") do
        context("called directly") do
            @fact Cₚ(1) => roughly(e)
            @fact Cₚ(0) => 1
            @fact_throws Cₚ(1,2)
        end

        context("called through a ValueSet") do
            vs1 = Ogre.ValueSet(1,2,3)
            vs2 = Ogre.ValueSet(1,2,3,4)

            @fact_throws Cₚ(vs1)
            @fact Cₚ(vs2) => exp(4)
        end
    end
end


