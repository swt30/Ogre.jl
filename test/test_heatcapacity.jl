include("header.jl")

facts("Heat capacity handling") do
    context("Construction") do
        Cₚ1 = Ogre.HeatCapacity(4200)
        Cₚ2 = Ogre.HeatCapacity(Ogre.WithTemp, T -> 4200)
        @fact Cₚ1(300) => Cₚ2(300)
        @fact_throws Ogre.HeatCapacity(Ogre.NoTemp, res.heat_capacity_function)
    end

    context("Values are as expected") do
        Cₚ = res.heatcap.exponential
        Cₚ2 = res.heatcap.exponential2d
        context("called directly") do
            @fact Cₚ(1) => roughly(e)
            @fact Cₚ(0) => 1
            @fact Cₚ2(0, 0) => 1
            @fact Cₚ2(1, 1) => roughly(e^2)
            @fact Cₚ(1,2) => roughly(e)
            @fact Cₚ(2,1) => roughly(e^2)
            @fact_throws Cₚ2(1)
        end

        context("called through a ValueSet") do
            @fact_throws Cₚ(res.value_set_no_temp)
            @fact Cₚ(res.value_set_full) => exp(4)
        end
    end
end

