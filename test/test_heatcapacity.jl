include("header.jl")

facts("Heat capacity handling") do
    capfunc(T) = exp(T)

    Cₚ = Ogre.HeatCapacity(Ogre.WithTemp, capfunc)
    Cₚ₂ = Ogre.HeatCapacity(capfunc)

    context("Construction") do
        @fact_throws Ogre.HeatCapacity(Ogre.NoTemp, capfunc)
    end

    context("Values are as expected") do
        context("called directly") do
            @fact Cₚ(1) => roughly(e)
            @fact Cₚ(0) => 1
            @fact_throws Cₚ(1,2)
        end

        context("called through a ValueSet") do
            @fact_throws Cₚ(res.value_set_no_temp)
            @fact Cₚ(res.value_set_full) => exp(4)
        end
    end
end


