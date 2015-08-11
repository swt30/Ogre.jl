include("header.jl")

facts("Heat capacity handling") do
    context("Values are as expected") do
        cₚ_exp = res.heatcap.exponential
        cₚ_exp2 = res.heatcap.exponential2d
        cₚ_const = res.heatcap.constant

        context("called directly") do
            # constant values are always equal
            @fact cₚ_const(0) --> cₚ_const(10)
            @fact cₚ_const(0, 10) --> 1

            # calling a single-input heat capacity
            @fact cₚ_exp(1) --> e^1
            @fact cₚ_exp(0) --> e^0
            @fact cₚ_exp(1, 2) --> e^2

            # calling a temp- and pressure-dependent heat capacity
            @fact cₚ_exp2(0, 1) --> e^0 + 2e^1
            @fact cₚ_exp2(1, 1) --> e^1 + 2e^1
            @fact cₚ_exp2(1, 3) --> e^1 + 2e^3
            @fact cₚ_exp2(2, 1) --> e^2 + 2e^1

            # needs both temperature and pressure
            @fact_throws MethodError cₚ_exp2(1)
        end

        context("called through a ValueSet") do
            # accessing temperature and pressure
            @fact_throws MethodError cₚ_exp(res.value_set_no_temp)
            @fact Ogre.temperature(res.value_set_full) --> 4
            @fact Ogre.pressure(res.value_set_full) --> 3

            # uses just temperature
            @fact cₚ_exp(res.value_set_full) --> e^4

            # uses temperature and pressure
            @fact cₚ_exp2(res.value_set_full) --> e^3 + 2e^4
        end
    end
end

